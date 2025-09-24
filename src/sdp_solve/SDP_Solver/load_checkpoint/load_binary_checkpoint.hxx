#pragma once

#include "sdpb_util/assert.hxx"
#include "sdpb_util/Verbosity.hxx"

#include <filesystem>
#include <boost/property_tree/json_parser.hpp>

// Helper functions to load binary data into DistMatrix etc.

inline void read_local_binary_matrix(El::AbstractDistMatrix<El::BigFloat> &m,
                                     std::ifstream &checkpoint_stream,
                                     std::vector<uint8_t> &buffer)
{
  int64_t local_height, local_width;
  checkpoint_stream.read(reinterpret_cast<char *>(&local_height),
                         sizeof(int64_t));
  checkpoint_stream.read(reinterpret_cast<char *>(&local_width),
                         sizeof(int64_t));

  ASSERT(checkpoint_stream.good(),
         "Corrupted binary checkpoint file.  For block with ", "global size (",
         m.Height(), ",", m.Width(), ") and local dimensions (",
         m.LocalHeight(), ",", m.LocalWidth(),
         "), error when reading height and width");

  ASSERT(local_height == m.LocalHeight() && local_width == m.LocalWidth(),
         "Incompatible binary checkpoint file.  For block with ",
         "global size (", m.Height(), ",", m.Width(),
         "), expected local dimensions (", m.LocalHeight(), ",",
         m.LocalWidth(), "), but found (", local_height, ",", local_width,
         ")");

  for(int64_t row = 0; row < local_height; ++row)
    for(int64_t column = 0; column < local_width; ++column)
      {
        auto &elt = m.Matrix().Ref(row, column);
        buffer.resize(elt.SerializedSize());
        checkpoint_stream.read(reinterpret_cast<char *>(buffer.data()),
                               std::streamsize(buffer.size()));
        ASSERT(checkpoint_stream.good(),
               "Corrupted binary checkpoint file. For block with ",
               "global size (", m.Height(), ",", m.Width(),
               ") and local dimensions (", m.LocalHeight(), ",",
               m.LocalWidth(), "), error when reading element (", row, ",",
               column, ")");
        elt.Deserialize(buffer.data());
      }
}

inline void read_local_binary_matrix(El::AbstractDistMatrix<El::BigFloat> &m,
                                     std::ifstream &checkpoint_stream)
{
  std::vector<uint8_t> buffer;
  read_local_binary_matrix(m, checkpoint_stream, buffer);
}

template <class TBlock_Matrix>
void read_local_binary_blocks(TBlock_Matrix &m,
                              std::ifstream &checkpoint_stream)
{
  std::vector<uint8_t> buffer;
  for(auto &block : m.blocks)
    read_local_binary_matrix(block, checkpoint_stream, buffer);
}

// Load binary checkpoint

template <class TSolver>
bool load_binary_checkpoint(const std::filesystem::path &checkpoint_directory,
                            const Verbosity &verbosity, TSolver &solver)
{
  namespace fs = std::filesystem;
  int64_t current_generation(-1), backup_generation(-1);
  if(El::mpi::Rank() == 0)
    {
      fs::path metadata(checkpoint_directory / "checkpoint.json");
      if(exists(metadata))
        {
          boost::property_tree::ptree tree;
          boost::property_tree::read_json(metadata.string(), tree);
          boost::optional<int64_t> current(
            tree.get_optional<int64_t>("current"));
          if(current)
            {
              current_generation = current.value();
            }
          else
            {
              RUNTIME_ERROR("Invalid or missing element 'current' in ",
                            metadata);
            }
          boost::optional<int64_t> backup(
            tree.get_optional<int64_t>("backup"));
          if(backup)
            {
              backup_generation = backup.value();
            }
        }
    }

  // We cast to El::byte.  Trying to use El::mpi::Broadcast() directly
  // with an int64_t leads to weird memory errors.  Perhaps it calls
  // the Broadcast() with int32_t instead of int64_t?
  El::mpi::Broadcast(reinterpret_cast<El::byte *>(&current_generation),
                     sizeof(current_generation) / sizeof(El::byte), 0,
                     El::mpi::COMM_WORLD);

  const auto get_checkpoint_path_for_rank
    = [&checkpoint_directory, &current_generation](const int mpi_rank) {
        const auto filename = current_generation == -1
                                ? "checkpoint." + std::to_string(mpi_rank)
                                : "checkpoint_"
                                    + std::to_string(current_generation) + "_"
                                    + std::to_string(mpi_rank);
        return checkpoint_directory / filename;
      };

  {
    // Check that all checkpoint files exist and there are no extra files.
    // all_files_exist = true if checkpoint files exist for all ranks=0..n-1
    // all_files_exist = false if no checkpoint files exist (for brevity, we check rank=0 only)
    // An error is thrown if there are extra files (we check for rank=n)
    El::byte all_files_exist = false;
    // Check only on rank=0 to avoid extra filesystem load.
    if(El::mpi::Rank() == 0)
      {
        // check rank=0
        const auto path_0 = get_checkpoint_path_for_rank(0);
        if(exists(path_0))
          {
            // check rank=1..n-1
            for(int rank = 1; rank < El::mpi::Size(); ++rank)
              {
                const auto curr_path = get_checkpoint_path_for_rank(rank);
                if(!exists(curr_path))
                  {
                    RUNTIME_ERROR("Cannot find checkpoint file for MPI rank=",
                                  rank, ": ", curr_path);
                  }
              }
            all_files_exist = true;

            // Check rank=n (should not exist)
            const auto out_of_bound_path
              = get_checkpoint_path_for_rank(El::mpi::Size());
            if(exists(out_of_bound_path))
              {
                RUNTIME_ERROR(
                  "Checkpoint should contain only files for MPI ranks from 0 "
                  "to ",
                  El::mpi::Size() - 1,
                  ", but found unexpected checkpoint file: ",
                  out_of_bound_path);
              }
          }

        if(!all_files_exist && current_generation != -1)
          {
            // Checkpoint should exist in this case
            RUNTIME_ERROR("Cannot find checkpoint file: ", path_0);
          }
      }
    El::mpi::Broadcast(all_files_exist, 0, El::mpi::COMM_WORLD);
    if(!all_files_exist && current_generation == -1)
      return false;
  }

  const fs::path checkpoint_file_path
    = get_checkpoint_path_for_rank(El::mpi::Rank());

  // Update generations
  if(current_generation != -1)
    {
      solver.current_generation = current_generation;
      // See note above about Broadcast()
      El::mpi::Broadcast(reinterpret_cast<El::byte *>(&backup_generation),
                         sizeof(current_generation) / sizeof(El::byte), 0,
                         El::mpi::COMM_WORLD);
    }
  else
    {
      current_generation = 0;
    }

  try
    {
      std::ifstream checkpoint_stream(checkpoint_file_path);
      ASSERT(checkpoint_stream.good(), "Cannot open ", checkpoint_file_path);
      if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
        {
          El::Output("Loading binary checkpoint from : ",
                     checkpoint_directory);
        }
      load_local_binary_data(solver, checkpoint_stream);
    }
  catch(std::exception &e)
    {
      RUNTIME_ERROR("Failed to read checkpoint from binary file ",
                    checkpoint_file_path, "\n  ", e.what());
    }
  solver.current_generation = current_generation;
  if(backup_generation != -1)
    {
      solver.backup_generation = backup_generation;
    }
  return true;
}
