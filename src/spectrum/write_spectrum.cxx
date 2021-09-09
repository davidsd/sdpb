#include "Zeros.hxx"
#include "../set_stream_precision.hxx"

#include <boost/filesystem.hpp>

void write_spectrum(const boost::filesystem::path &output_path,
                    const std::vector<Zeros> &zeros_blocks)
{
  // Synchronize zeros
  // This is waaaay more work than it should be.
  std::vector<int32_t> sizes(zeros_blocks.size(), 0);
  for(size_t index(0); index != zeros_blocks.size(); ++index)
    {
      sizes.at(index) = zeros_blocks.at(index).zeros.size();
    }
  El::mpi::AllReduce(sizes.data(), sizes.size(), El::mpi::MAX,
                     El::mpi::COMM_WORLD);
  const int64_t num_elements(std::accumulate(sizes.begin(), sizes.end(), 0));
  std::vector<El::BigFloat> zero_flattened;
  zero_flattened.reserve(num_elements);
  for(size_t block_index(0); block_index != zeros_blocks.size(); ++block_index)
    {
      for(int32_t zero_index(0); zero_index != sizes[block_index];
          ++zero_index)
        {
          if(zeros_blocks.at(block_index).zeros.empty())
            {
              zero_flattened.emplace_back(0);
            }
          else
            {
              zero_flattened.emplace_back(
                zeros_blocks.at(block_index).zeros.at(zero_index));
            }
        }
    }
  El::mpi::Reduce(zero_flattened.data(), zero_flattened.size(), 0,
                  El::mpi::COMM_WORLD);

  // Write the file
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
      auto iterator(zero_flattened.begin());
      for(size_t block_index(0); block_index != sizes.size(); ++block_index)
        {
          if(block_index != 0)
            {
              outfile << ",";
            }
          outfile << "\n  [";
          for(int64_t zero_index(0); zero_index != sizes.at(block_index);
              ++zero_index)
            {
              if(zero_index != 0)
                {
                  outfile << ",";
                }
              outfile << "\n    \"" << *iterator << "\"";
              ++iterator;
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
}
