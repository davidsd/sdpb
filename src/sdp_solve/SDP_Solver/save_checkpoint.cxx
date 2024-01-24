#include "../SDP_Solver.hxx"
#include "sdpb_util/assert.hxx"

#include <filesystem>
#include <boost/property_tree/json_parser.hpp>

namespace fs = std::filesystem;

// We use binary checkpointing because writing text does not write all
// of the necessary digits.  The GMP library sets it to one less than
// required for round-tripping.
template <typename T>
void write_local_blocks(const T &t, std::ofstream &checkpoint_stream)
{
  El::BigFloat zero(0);
  const size_t serialized_size(zero.SerializedSize());
  std::vector<uint8_t> local_array(serialized_size);

  for(auto &block : t.blocks)
    {
      int64_t local_height(block.LocalHeight()),
        local_width(block.LocalWidth());
      checkpoint_stream.write(reinterpret_cast<char *>(&local_height),
                              sizeof(int64_t));
      checkpoint_stream.write(reinterpret_cast<char *>(&local_width),
                              sizeof(int64_t));
      for(int64_t row = 0; row < local_height; ++row)
        for(int64_t column = 0; column < local_width; ++column)
          {
            block.GetLocal(row, column).Serialize(local_array.data());
            checkpoint_stream.write(
              reinterpret_cast<char *>(local_array.data()),
              std::streamsize(local_array.size()));
          }
    }
}

void SDP_Solver::save_checkpoint(
  const fs::path &checkpoint_directory, const Verbosity &verbosity,
  const boost::property_tree::ptree &parameter_properties)
{
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
  if(backup_generation)
    {
      remove(checkpoint_directory
             / ("checkpoint_" + std::to_string(backup_generation.value()) + "_"
                + std::to_string(El::mpi::Rank())));
    }
  backup_generation = current_generation;
  current_generation += 1;
  fs::path checkpoint_filename(checkpoint_directory
    / ("checkpoint_" + std::to_string(current_generation) + "_"
       + std::to_string(El::mpi::Rank())));

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
      // TODO: Write and read precision, num of mpi procs, and procs_per_node.
      write_local_blocks(x, checkpoint_stream);
      write_local_blocks(X, checkpoint_stream);
      write_local_blocks(y, checkpoint_stream);
      write_local_blocks(Y, checkpoint_stream);
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
      metadata << "{\n    \"current\": " << current_generation << ",\n"
               << "    \"backup\": " << backup_generation.value() << ",\n"
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
