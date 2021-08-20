#include "../set_stream_precision.hxx"

#include <El.hpp>

#include <boost/filesystem.hpp>

void write_spectrum(const boost::filesystem::path &output_path,
                    const std::vector<std::vector<El::BigFloat>> &zeros)
{
  boost::filesystem::ofstream outfile(output_path);
  if(!outfile.good())
    {
      throw std::runtime_error("Problem when opening output file: '"
                               + output_path.string() + "'");
    }
  set_stream_precision(outfile);
  outfile << "[";
  for(auto block(zeros.begin()); block != zeros.end(); ++block)
    {
      if(block != zeros.begin())
        {
          outfile << ",";
        }
      outfile << "\n  [";
      for(auto zero(block->begin()); zero != block->end(); ++zero)
        {
          if(zero != block->begin())
            {
              outfile << ",";
            }
          outfile << "\n    \"" << *zero << "\"";
        }
      outfile << "\n  ]";
    }
  outfile << "\n]\n";
  if(!outfile.good())
    {
      throw std::runtime_error("Problem when writing to output file: '"
                               + output_path.string() + "'");
    }
}
