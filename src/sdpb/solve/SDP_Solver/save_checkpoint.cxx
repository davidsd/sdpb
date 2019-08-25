#include "../SDP_Solver.hxx"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

// We use binary checkpointing because writing text does not write all
// of the necessary digits.  The GMP library sets it to one less than
// required for round-tripping.
template <typename T>
void write_local_blocks(const T &t,
                        boost::filesystem::ofstream &checkpoint_stream)
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
  const boost::filesystem::path &checkpoint_directory,
  const Verbosity &verbosity) const
{
  boost::filesystem::path checkpoint_filename(
    checkpoint_directory / ("checkpoint." + std::to_string(El::mpi::Rank())));

  if(!exists(checkpoint_directory))
    {
      create_directories(checkpoint_directory);
    }
  else if(!is_directory(checkpoint_directory))
    {
      throw std::runtime_error("Checkpoint directory '"
                               + checkpoint_directory.string()
                               + "'already exists, but is not a directory");
    }
  if(exists(checkpoint_directory / checkpoint_filename))
    {
      if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
        {
          std::cout << "Backing up checkpoint\n";
        }
      boost::filesystem::path backup_filename(checkpoint_filename.string()
                                              + ".bk");
      remove(backup_filename);
      rename(checkpoint_filename, backup_filename);
    }

  const size_t max_retries(10);
  bool wrote_successfully(false);
  for(size_t attempt=0; attempt<max_retries && !wrote_successfully; ++attempt)
    {
      boost::filesystem::ofstream checkpoint_stream(checkpoint_filename);
      if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
        {
          std::cout << "Saving checkpoint to    : " << checkpoint_directory
                    << '\n';
        }
      write_local_blocks(x, checkpoint_stream);
      write_local_blocks(X, checkpoint_stream);
      write_local_blocks(y, checkpoint_stream);
      write_local_blocks(Y, checkpoint_stream);
      wrote_successfully=checkpoint_stream.good();
      if(!wrote_successfully)
        {
          if(attempt+1<max_retries)
            {
              std::stringstream ss;
              ss << "Error writing checkpoint file '"
                 << checkpoint_filename << "'.  Retrying "
                 << (attempt+2) << "/" << max_retries << "\n";
              std::cerr << ss.str() << std::flush;
            }
          else
            {
              std::stringstream ss;
              ss << "Error writing checkpoint file '"
                 << checkpoint_filename << "'.  Exceeded max retries.\n";
              throw std::runtime_error(ss.str());
            }
        }
    }
}
