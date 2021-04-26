#include <El.hpp>

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
