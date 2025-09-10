#pragma once

#include "sdpb_util/assert.hxx"
#include "sdpb_util/Verbosity.hxx"

#include <filesystem>
#include <boost/property_tree/json_parser.hpp>

// We use binary checkpointing because writing text does not write all
// the necessary digits.  The GMP library sets it to one less than
// required for round-tripping.
inline void write_local_matrix(const El::AbstractDistMatrix<El::BigFloat> &m,
                               std::ofstream &checkpoint_stream,
                               std::vector<uint8_t> &buffer)
{
  int64_t local_height(m.LocalHeight()), local_width(m.LocalWidth());
  checkpoint_stream.write(reinterpret_cast<char *>(&local_height),
                          sizeof(int64_t));
  checkpoint_stream.write(reinterpret_cast<char *>(&local_width),
                          sizeof(int64_t));
  for(int64_t row = 0; row < local_height; ++row)
    for(int64_t column = 0; column < local_width; ++column)
      {
        const auto &elt = m.GetLocalCRef(row, column);
        buffer.resize(elt.SerializedSize());
        elt.Serialize(buffer.data());
        checkpoint_stream.write(reinterpret_cast<char *>(buffer.data()),
                                std::streamsize(buffer.size()));
      }
}

inline void write_local_matrix(const El::AbstractDistMatrix<El::BigFloat> &m,
                               std::ofstream &checkpoint_stream)
{
  std::vector<uint8_t> buffer;
  write_local_matrix(m, checkpoint_stream, buffer);
}
template <class TBlock_Matrix>
void write_local_blocks(const TBlock_Matrix &m,
                        std::ofstream &checkpoint_stream)
{
  std::vector<uint8_t> buffer;
  for(auto &block : m.blocks)
    {
      write_local_matrix(block, checkpoint_stream, buffer);
    }
}

template <class TSolver>
void save_binary_checkpoint(
  const std::filesystem::path &checkpoint_directory,
  const Verbosity &verbosity,
  const boost::property_tree::ptree &parameter_properties, TSolver &solver)
{
  namespace fs = std::filesystem;
  if(checkpoint_directory.empty())
    {
      return;
    }
  if(!exists(checkpoint_directory))
    {
      create_directories(checkpoint_directory);
    }
  else if(!is_directory(checkpoint_directory))
    {
      RUNTIME_ERROR(
        "Checkpoint directory already exists, but is not a directory: ",
        checkpoint_directory);
    }
  if(solver.backup_generation)
    {
      remove(checkpoint_directory
             / ("checkpoint_"
                + std::to_string(solver.backup_generation.value()) + "_"
                + std::to_string(El::mpi::Rank())));
    }
  solver.backup_generation = solver.current_generation;
  solver.current_generation += 1;
  fs::path checkpoint_filename(checkpoint_directory
                               / ("checkpoint_"
                                  + std::to_string(solver.current_generation)
                                  + "_" + std::to_string(El::mpi::Rank())));

  const size_t max_retries(10);
  bool wrote_successfully(false);
  for(size_t attempt = 0; attempt < max_retries && !wrote_successfully;
      ++attempt)
    {
      std::ofstream checkpoint_stream(checkpoint_filename);
      if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
        {
          std::cout << "Saving checkpoint to    : " << checkpoint_directory
                    << '\n';
        }
      write_binary_checkpoint_data(checkpoint_stream, solver);
      wrote_successfully = checkpoint_stream.good();
      if(!wrote_successfully)
        {
          if(attempt + 1 < max_retries)
            {
              std::stringstream ss;
              ss << "Error writing checkpoint file '" << checkpoint_filename
                 << "'.  Retrying " << (attempt + 2) << "/" << max_retries
                 << "\n";
              std::cerr << ss.str() << std::flush;
            }
          else
            {
              RUNTIME_ERROR("Error writing checkpoint file ",
                            checkpoint_filename, ":  Exceeded max retries.");
            }
        }
    }
  if(El::mpi::Rank() == 0)
    {
      std::ofstream metadata(checkpoint_directory / "checkpoint_new.json");
      metadata << "{\n    \"current\": " << solver.current_generation << ",\n"
               << "    \"backup\": " << solver.backup_generation.value()
               << ",\n"
               << "    \"version\": \"" << SDPB_VERSION_STRING
               << "\",\n    \"options\": \n";

      boost::property_tree::write_json(metadata, parameter_properties);
      metadata << "}\n";
    }
  El::mpi::Barrier(El::mpi::COMM_WORLD);
  if(El::mpi::Rank() == 0)
    {
      rename(checkpoint_directory / "checkpoint_new.json",
             checkpoint_directory / "checkpoint.json");
    }
}
