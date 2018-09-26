#include "../SDP_Solver.hxx"

#include <boost/archive/text_iarchive.hpp>

template <typename T>
void read_local_blocks(T &t,
                       boost::filesystem::ifstream &checkpoint_stream)
{
  for(auto &block : t.blocks)
    {
      int64_t local_height, local_width;
      checkpoint_stream >> local_height >> local_width;
      if(local_height != block.LocalHeight()
         || local_width != block.LocalWidth())
        {
          std::stringstream ss;
          ss << "Incompatible checkpoint file.  Expected dimensions ("
             << block.LocalHeight() << "," << block.LocalWidth()
             << "), but found (" << local_height << "," << local_width << ")";

          throw std::runtime_error(ss.str());
        }

      for(int64_t row = 0; row < local_height; ++row)
        for(int64_t column = 0; column < local_width; ++column)
          {
            El::BigFloat input;
            checkpoint_stream >> input;
            block.SetLocal(row, column, input);
          }
    }
}

void SDP_Solver::load_checkpoint(
  const boost::filesystem::path &checkpoint_directory)
{
  boost::filesystem::path checkpoint_filename(
    checkpoint_directory / ("checkpoint." + std::to_string(El::mpi::Rank())));

  if(!exists(checkpoint_filename))
    {
      return;
    }

  boost::filesystem::ifstream checkpoint_stream(checkpoint_filename);
  if(El::mpi::Rank() == 0)
    {
      std::cout << "Loading checkpoint from : " << checkpoint_directory
                << '\n';
    }
  read_local_blocks(x, checkpoint_stream);
  read_local_blocks(X, checkpoint_stream);
  read_local_blocks(y, checkpoint_stream);
  read_local_blocks(Y, checkpoint_stream);
}
