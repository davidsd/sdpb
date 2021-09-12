#include "../Zeros.hxx"

#include <boost/filesystem.hpp>

void write_file(const boost::filesystem::path &output_path,
                const std::vector<Zeros> &zeros_blocks);

void write_spectrum(const boost::filesystem::path &output_path,
                    const size_t &num_blocks,
                    const std::vector<Zeros> &zeros_blocks)
{
  // Synchronize zeros
  // This is waaaay more work than it should be.
  const size_t rank(El::mpi::Rank()),
    num_procs(El::mpi::Size(El::mpi::COMM_WORLD));

  if(num_procs == 1)
    {
      write_file(output_path, zeros_blocks);
    }
  else
    {
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
      // El::mpi::Reduce(m.Buffer(), m.Height()*m.Width(), 0,
      // El::mpi::COMM_WORLD);

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
      write_file(output_path, zeros_all_blocks);
    }
}
