#include "sdpb_util/assert.hxx"
#include "spectrum/Zeros.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <filesystem>

namespace fs = std::filesystem;

void write_file(const fs::path &output_path,
                const std::vector<Zeros> &zeros_blocks)
{
  if(El::mpi::Rank() == 0)
    {
      if(output_path.has_parent_path())
        fs::create_directories(output_path.parent_path());
      std::ofstream outfile(output_path);
      ASSERT(outfile.good(),
             "Problem when opening output file: ", output_path);
      set_stream_precision(outfile);
      outfile << "[";
      for(auto zeros_iterator(zeros_blocks.begin());
          zeros_iterator != zeros_blocks.end(); ++zeros_iterator)
        {
          if(zeros_iterator != zeros_blocks.begin())
            {
              outfile << ",";
            }
          auto block_path = zeros_iterator->block_path.string();
          ASSERT(!block_path.empty(), "Empty path for block_",
                 std::distance(zeros_blocks.begin(), zeros_iterator));
          outfile << "\n  {\n    \"block_path\": \"" << block_path << "\",";
          outfile << "\n    \"zeros\":\n      [";
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
      ASSERT(outfile.good(),
             "Problem when writing to output file: ", output_path);
    }
}
