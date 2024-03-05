#include "sdp_solve/SDP_Solver.hxx"
#include "sdpb_util/assert.hxx"

#include <filesystem>
#include <boost/property_tree/json_parser.hpp>

namespace fs = std::filesystem;

template <typename T>
void read_local_binary_blocks(T &t, std::ifstream &checkpoint_stream)
{
  El::BigFloat zero(0);
  const size_t serialized_size(zero.SerializedSize());
  std::vector<uint8_t> local_array(serialized_size);

  for(auto &block : t.blocks)
    {
      int64_t local_height, local_width;
      checkpoint_stream.read(reinterpret_cast<char *>(&local_height),
                             sizeof(int64_t));
      checkpoint_stream.read(reinterpret_cast<char *>(&local_width),
                             sizeof(int64_t));

      ASSERT(checkpoint_stream.good(),
             "Corrupted binary checkpoint file.  For block with ",
             "global size (", block.Height(), ",", block.Width(),
             ") and local dimensions (", block.LocalHeight(), ",",
             block.LocalWidth(), "), error when reading height and width");

      ASSERT(local_height == block.LocalHeight()
               && local_width == block.LocalWidth(),
             "Incompatible binary checkpoint file.  For block with ",
             "global size (", block.Height(), ",", block.Width(),
             "), expected local dimensions (", block.LocalHeight(), ",",
             block.LocalWidth(), "), but found (", local_height, ",",
             local_width, ")");

      for(int64_t row = 0; row < local_height; ++row)
        for(int64_t column = 0; column < local_width; ++column)
          {
            El::BigFloat input;
            checkpoint_stream.read(
              reinterpret_cast<char *>(local_array.data()),
              std::streamsize(local_array.size()));
            ASSERT(checkpoint_stream.good(),
                   "Corrupted binary checkpoint file. For block with ",
                   "global size (", block.Height(), ",", block.Width(),
                   ") and local dimensions (", block.LocalHeight(), ",",
                   block.LocalWidth(), "), error when reading element (", row,
                   ",", column, ")");
            input.Deserialize(local_array.data());

            block.SetLocal(row, column, input);
          }
    }
}

bool load_binary_checkpoint(const fs::path &checkpoint_directory,
                            const Verbosity &verbosity, SDP_Solver &solver)
{
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
  read_local_binary_blocks(solver.x, checkpoint_stream);
  read_local_binary_blocks(solver.X, checkpoint_stream);
  read_local_binary_blocks(solver.y, checkpoint_stream);
  read_local_binary_blocks(solver.Y, checkpoint_stream);
  solver.current_generation = current_generation;
  if(backup_generation != -1)
    {
      solver.backup_generation = backup_generation;
    }
  return true;
}
