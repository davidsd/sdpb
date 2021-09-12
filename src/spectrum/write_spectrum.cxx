#include "Zeros.hxx"
#include "../set_stream_precision.hxx"

#include <boost/filesystem.hpp>

void write_spectrum(const boost::filesystem::path &output_path,
                    const size_t &num_blocks,
                    const std::vector<Zeros> &zeros_blocks)
{
  // Synchronize zeros
  // This is waaaay more work than it should be.
  const size_t rank(El::mpi::Rank()),
    num_procs(El::mpi::Size(El::mpi::COMM_WORLD));

  std::vector<size_t> zero_sizes(num_blocks, 0);
  size_t block_index(rank);
  for(auto &block : zeros_blocks)
    {
      zero_sizes.at(block_index) = block.zeros.size();
      block_index += num_procs;
    }
  El::mpi::AllReduce(zero_sizes.data(), zero_sizes.size(), El::mpi::SUM,
                     El::mpi::COMM_WORLD);

  // El::Matrix<El::BigFloat> m(10,10);
  // El::Fill(m,El::BigFloat(1));
  // El::mpi::Reduce(m.Buffer(), m.Height()*m.Width(), 0, El::mpi::COMM_WORLD);

  std::vector<Zeros> zeros_all_blocks(num_blocks);
  for(size_t block_index(0); block_index != zeros_all_blocks.size();
      ++block_index)
    {
      if(block_index % num_procs == rank)
        {
          const size_t local_index(block_index / num_procs);
          zeros_all_blocks[block_index].zeros
            = zeros_blocks.at(local_index).zeros;
        }
      else
        {
          zeros_all_blocks[block_index].zeros.resize(
            zero_sizes.at(block_index), El::BigFloat(0));
        }
      El::mpi::Reduce(zeros_all_blocks[block_index].zeros.data(),
                      zeros_all_blocks[block_index].zeros.size(), 0,
                      El::mpi::COMM_WORLD);
    }

  // Write the file
  if(El::mpi::Rank() == 0)
    {
      // El::Print(m,"m");
      // std::cout << "\n";
      boost::filesystem::ofstream outfile(output_path);
      if(!outfile.good())
        {
          throw std::runtime_error("Problem when opening output file: '"
                                   + output_path.string() + "'");
        }
      set_stream_precision(outfile);
      outfile << "[";
      for(auto zeros_iterator(zeros_all_blocks.begin());
          zeros_iterator != zeros_all_blocks.end(); ++zeros_iterator)
        {
          if(zeros_iterator != zeros_all_blocks.begin())
            {
              outfile << ",";
            }
          outfile << "\n  [";
          for(size_t zero_index(0); zero_index != zeros_iterator->zeros.size();
              ++zero_index)
            {
              if(zero_index != 0)
                {
                  outfile << ",";
                }
              outfile << "\n    \"" << zeros_iterator->zeros.at(zero_index)
                      << "\"";
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
