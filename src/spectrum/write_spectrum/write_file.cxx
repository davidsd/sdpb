#include "../Zeros.hxx"
#include "../../set_stream_precision.hxx"

#include <boost/filesystem.hpp>

void write_file(const boost::filesystem::path &output_path,
                const std::vector<Zeros> &zeros_blocks)
{
  if(El::mpi::Rank() == 0)
    {
      boost::filesystem::ofstream outfile(output_path);
      if(!outfile.good())
        {
          throw std::runtime_error("Problem when opening output file: '"
                                   + output_path.string() + "'");
        }
      set_stream_precision(outfile);
      outfile << "[";
      for(auto zeros_iterator(zeros_blocks.begin());
          zeros_iterator != zeros_blocks.end(); ++zeros_iterator)
        {
          if(zeros_iterator != zeros_blocks.begin())
            {
              outfile << ",";
            }
          outfile << "\n  {\n    \"zeros\":\n      [";
          for(size_t zero_index(0); zero_index != zeros_iterator->zeros.size();
              ++zero_index)
            {
              if(zero_index != 0)
                {
                  outfile << ",";
                }
              outfile << "\n        {\n          \"zero\": \""
                      << zeros_iterator->zeros.at(zero_index).zero << "\""
                      << ",\n          \"lambda\":\n            [\n";
              for(int64_t row(0);
                  row < zeros_iterator->zeros.at(zero_index).lambda.Height();
                  ++row)
                {
                  if(row != 0)
                    {
                      outfile << ",\n";
                    }
                  outfile << "              \""
                          << zeros_iterator->zeros.at(zero_index).lambda(row, 0)
                          << "\"";
                }
              outfile << "\n            ]\n        }";
            }
          outfile << "\n      ],\n"
                  << "    \"error\": \"" << zeros_iterator->error << "\"\n"
                  << "  }";
        }
      outfile << "\n]\n";
      if(!outfile.good())
        {
          throw std::runtime_error("Problem when writing to output file: '"
                                   + output_path.string() + "'");
        }
    }
}
