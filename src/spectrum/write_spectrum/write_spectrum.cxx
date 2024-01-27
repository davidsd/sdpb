#include "sdpb_util/assert.hxx"
#include "spectrum/Zeros.hxx"

#include <filesystem>

namespace fs = std::filesystem;

void write_file(const fs::path &output_path,
                const std::vector<Zeros> &zeros_blocks);

void write_spectrum(const fs::path &output_path, const size_t &num_blocks,
                    const std::vector<Zeros> &zeros_blocks,
                    const std::vector<size_t> &block_indices)
{
  if(El::mpi::Size() == 1)
    {
      write_file(output_path, zeros_blocks);
    }
  else
    {
      // Synchronize zeros
      // This is waaaay more work than it should be.
      std::vector<size_t> zero_sizes(num_blocks, 0),
        lambda_sizes(num_blocks, 0);

      ASSERT(block_indices.size() == zeros_blocks.size(),
             "block_indices.size()=", block_indices.size(),
             " and zeros_blocks.size()=", zeros_blocks.size(),
             " should be equal");
      for(size_t local_index = 0; local_index < block_indices.size();
          ++local_index)
        {
          const auto &block = zeros_blocks.at(local_index);
          size_t block_index = block_indices.at(local_index);
          zero_sizes.at(block_index) = block.zeros.size();
          if(!block.zeros.empty())
            {
              lambda_sizes.at(block_index)
                = block.zeros.front().lambda.Height();
            }
        }
      El::mpi::AllReduce(zero_sizes.data(), zero_sizes.size(), El::mpi::SUM,
                         El::mpi::COMM_WORLD);
      El::mpi::AllReduce(lambda_sizes.data(), lambda_sizes.size(),
                         El::mpi::SUM, El::mpi::COMM_WORLD);

      std::vector<Zeros> zeros_all_blocks(num_blocks);
      for(size_t block_index(0); block_index != num_blocks; ++block_index)
        {
          // TODO this search is not optimal, O(N) for each block_index
          auto local_index_it = std::find(block_indices.begin(),
                                          block_indices.end(), block_index);
          if(local_index_it != block_indices.end())
            {
              size_t local_index = local_index_it - block_indices.begin();
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
