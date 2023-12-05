#include "spectrum/Zeros.hxx"

#include <filesystem>

namespace fs = std::filesystem;

void write_file(const fs::path &output_path,
                const std::vector<Zeros> &zeros_blocks);

void write_spectrum(const fs::path &output_path, const size_t &num_blocks,
                    const std::vector<Zeros> &zeros_blocks)
{
  const size_t rank(El::mpi::Rank()),
    num_procs(El::mpi::Size(El::mpi::COMM_WORLD));

  if(num_procs == 1)
    {
      write_file(output_path, zeros_blocks);
    }
  else
    {
      // Synchronize zeros
      // This is waaaay more work than it should be.
      std::vector<size_t> zero_sizes(num_blocks, 0),
        lambda_sizes(num_blocks, 0);
      size_t block_index(rank);
      for(auto &block : zeros_blocks)
        {
          zero_sizes.at(block_index) = block.zeros.size();
          if(!block.zeros.empty())
            {
              lambda_sizes.at(block_index)
                = block.zeros.front().lambda.Height();
            }
          block_index += num_procs;
        }
      El::mpi::AllReduce(zero_sizes.data(), zero_sizes.size(), El::mpi::SUM,
                         El::mpi::COMM_WORLD);
      El::mpi::AllReduce(lambda_sizes.data(), lambda_sizes.size(),
                         El::mpi::SUM, El::mpi::COMM_WORLD);

      std::vector<Zeros> zeros_all_blocks(num_blocks);
      for(size_t block_index(0); block_index != zeros_all_blocks.size();
          ++block_index)
        {
          if(block_index % num_procs == rank)
            {
              const size_t local_index(block_index / num_procs);
              zeros_all_blocks[block_index] = zeros_blocks.at(local_index);
            }
          else
            {
              zeros_all_blocks[block_index].zeros.resize(
                zero_sizes.at(block_index), El::BigFloat(0));
              for(auto &zero : zeros_all_blocks[block_index].zeros)
                {
                  zero.lambda.Resize(lambda_sizes.at(block_index), 1);
                  El::Zero(zero.lambda);
                }
              zeros_all_blocks[block_index].error = 0;
            }

          for(auto &zero : zeros_all_blocks[block_index].zeros)
            {
              // TODO: This is painfully slow with lots of communication.
              zero.zero = El::mpi::Reduce(zero.zero, El::mpi::SUM, 0,
                                          El::mpi::COMM_WORLD);
              El::mpi::Reduce(zero.lambda.Buffer(),
                              zero.lambda.Height() * zero.lambda.Width(),
                              El::mpi::SUM, 0, El::mpi::COMM_WORLD);
            }
          zeros_all_blocks[block_index].error
            = El::mpi::Reduce(zeros_all_blocks[block_index].error,
                              El::mpi::SUM, 0, El::mpi::COMM_WORLD);
        }
      write_file(output_path, zeros_all_blocks);
    }
}
