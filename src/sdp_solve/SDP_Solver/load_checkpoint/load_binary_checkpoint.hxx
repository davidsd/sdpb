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
  fs::path checkpoint_filename;
  if(current_generation != -1)
    {
      solver.current_generation = current_generation;
      checkpoint_filename
        = checkpoint_directory
          / ("checkpoint_" + std::to_string(current_generation) + "_"
             + std::to_string(El::mpi::Rank()));
      ASSERT(exists(checkpoint_filename),
             "Missing checkpoint file: ", checkpoint_filename);
      // See note above about Broadcast()
      El::mpi::Broadcast(reinterpret_cast<El::byte *>(&backup_generation),
                         sizeof(current_generation) / sizeof(El::byte), 0,
                         El::mpi::COMM_WORLD);
    }
  else
    {
      checkpoint_filename
        = checkpoint_directory
          / ("checkpoint." + std::to_string(El::mpi::Rank()));
      if(!exists(checkpoint_filename))
        {
          return false;
        }
      current_generation = 0;
    }

  std::ifstream checkpoint_stream(checkpoint_filename);
  if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      std::cout << "Loading binary checkpoint from : " << checkpoint_directory
                << '\n';
    }
  load_local_binary_data(solver, checkpoint_stream);
  solver.current_generation = current_generation;
  if(backup_generation != -1)
    {
      solver.backup_generation = backup_generation;
    }
  return true;
}
