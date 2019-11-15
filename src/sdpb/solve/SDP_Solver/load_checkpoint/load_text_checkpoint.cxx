#include "../../SDP_Solver.hxx"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

void read_text_block(El::DistMatrix<El::BigFloat> &block,
                     const boost::filesystem::path &block_path)
{
  boost::filesystem::ifstream block_stream(block_path);
  if(!block_stream)
    {
      throw std::runtime_error("Unable to open checkpoint file: '"
                               + block_path.string() + "'");
    }
  int64_t file_height, file_width;
  block_stream >> file_height >> file_width;
  if(!block_stream.good())
    {
      throw std::runtime_error("Corrupted header in file: "
                               + block_path.string());
    }
  if(file_height != block.Height() || file_width != block.Width())
    {
      std::stringstream ss;
      ss << "Incompatible checkpoint file: '" << block_path.string()
         << "'.  Expected dimensions (" << block.Height() << ","
         << block.Width() << "), but found (" << file_height << ","
         << file_width << ")";
      throw std::runtime_error(ss.str());
    }

  std::string element;
  for(int64_t row = 0; row < file_height; ++row)
    for(int64_t column = 0; column < file_width; ++column)
      {
        block_stream >> element;
        if(block.IsLocal(row, column))
          {
            block.SetLocal(block.LocalRow(row), block.LocalCol(column),
                           El::BigFloat(element));
          }
      }
  if(!block_stream.good())
    {
      throw std::runtime_error("Corrupted data in file: "
                               + block_path.string());
    }
}

void read_text_block(El::DistMatrix<El::BigFloat> &block,
                     const boost::filesystem::path &checkpoint_directory,
                     const std::string &prefix, const size_t &block_index)
{
  read_text_block(block, checkpoint_directory
                           / (prefix + std::to_string(block_index) + ".txt"));
}

bool load_text_checkpoint(const boost::filesystem::path &checkpoint_directory,
                          const std::vector<size_t> &block_indices,
                          const Verbosity &verbosity, SDP_Solver &solver)
{
  if(!exists(checkpoint_directory / "x_0.txt"))
    {
      return false;
    }

  if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      std::cout << "Loading text checkpoint from : " << checkpoint_directory
                << '\n';
    }

  for(size_t block = 0; block != block_indices.size(); ++block)
    {
      size_t block_index(block_indices.at(block));
      read_text_block(solver.x.blocks.at(block), checkpoint_directory, "x_",
                      block_index);
      read_text_block(solver.y.blocks.at(block),
                      checkpoint_directory / "y.txt");

      for(size_t psd_block(0); psd_block < 2; ++psd_block)
        {
          // Constant constraints have empty odd parity blocks, so we do not
          // need to load them.
          if(solver.X.blocks.at(2 * block + psd_block).Height() != 0)
            {
              const size_t psd_index(2 * block_index + psd_block);
              read_text_block(solver.X.blocks.at(2 * block + psd_block),
                              checkpoint_directory, "X_matrix_", psd_index);
              read_text_block(solver.Y.blocks.at(2 * block + psd_block),
                              checkpoint_directory, "Y_matrix_", psd_index);
            }
        }
    }
  return true;
}
