#include "../SDP_Solver.hxx"
#include "../../Timers.hxx"
#include "../../../set_stream_precision.hxx"

#include <boost/filesystem.hpp>

template <typename T>
void write_local_blocks(const T &t,
                        boost::filesystem::ofstream &checkpoint_stream)
{
  for(auto &block : t.blocks)
    {
      checkpoint_stream << block.LocalHeight() << " " << block.LocalWidth()
                        << "\n";
      for(int64_t row = 0; row < block.LocalHeight(); ++row)
        for(int64_t column = 0; column < block.LocalWidth(); ++column)
          {
            checkpoint_stream << block.GetLocal(row, column) << "\n";
          }
    }
}

void SDP_Solver::save_checkpoint(
  const boost::filesystem::path &checkpoint_directory)
{
  boost::filesystem::path checkpoint_filename(
    checkpoint_directory / ("checkpoint." + std::to_string(El::mpi::Rank())));

  if(!exists(checkpoint_directory))
    {
      create_directory(checkpoint_directory);
    }
  else if(!is_directory(checkpoint_directory))
    {
      throw std::runtime_error("Checkpoint directory '"
                               + checkpoint_directory.string()
                               + "'already exists, but is not a directory");
    }
  if(exists(checkpoint_directory / checkpoint_filename))
    {
      if(El::mpi::Rank() == 0)
        {
          std::cout << "Backing up checkpoint\n";
        }
      boost::filesystem::path backup_filename(checkpoint_filename.string() + ".bk");
      remove(backup_filename);
      rename(checkpoint_filename, backup_filename);
    }
  boost::filesystem::ofstream checkpoint_stream(checkpoint_filename);
  if(El::mpi::Rank() == 0)
    {
      std::cout << "Saving checkpoint to    : " << checkpoint_directory
                << '\n';
    }
  set_stream_precision(checkpoint_stream);
  write_local_blocks(x, checkpoint_stream);
  write_local_blocks(X, checkpoint_stream);
  write_local_blocks(y, checkpoint_stream);
  write_local_blocks(Y, checkpoint_stream);
}
