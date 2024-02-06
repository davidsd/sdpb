#include "sdpb_util/assert.hxx"
#include "spectrum/Zeros.hxx"

#include <filesystem>

namespace fs = std::filesystem;

void write_file(const fs::path &output_path,
                const std::vector<Zeros> &zeros_blocks);

namespace
{
  // Send Zeros from some rank to rank=0
  void synchronize_zeros_block(Zeros &zeros_block, const int from)
  {
    const int to = 0;

    const int rank = El::mpi::Rank();
    if(rank != to && rank != from)
      return;
    if(to == from)
      return;

    // zeros
    {
      auto &zeros = zeros_block.zeros;

      // Synchronize zeros.size()
      if(rank == from)
        {
          El::mpi::Send<size_t>(zeros.size(), to, El::mpi::COMM_WORLD);
        }
      if(rank == to)
        {
          const auto zeros_size
            = El::mpi::Recv<size_t>(from, El::mpi::COMM_WORLD);
          zeros.resize(zeros_size, El::BigFloat(0));
        }

      // Synchronize zeros vector
      for(auto &zero : zeros)
        {
          // zero
          {
            if(rank == from)
              El::mpi::Send(zero.zero, to, El::mpi::COMM_WORLD);
            if(rank == to)
              zero.zero
                = El::mpi::Recv<El::BigFloat>(from, El::mpi::COMM_WORLD);
          }

          // lambda
          {
            auto &lambda = zero.lambda;
            if(rank == from)
              {
                El::mpi::Send(lambda.Height(), to, El::mpi::COMM_WORLD);
                El::mpi::Send(lambda.Width(), to, El::mpi::COMM_WORLD);
                El::Send(lambda, El::mpi::COMM_WORLD, to);
              }
            if(rank == to)
              {
                int height = El::mpi::Recv<int>(from, El::mpi::COMM_WORLD);
                int width = El::mpi::Recv<int>(from, El::mpi::COMM_WORLD);
                lambda.Resize(height, width);
                El::Recv(lambda, El::mpi::COMM_WORLD, from);
              }
          }
        }
    }

    // error
    {
      if(rank == from)
        El::mpi::Send(zeros_block.error, to, El::mpi::COMM_WORLD);
      if(rank == to)
        zeros_block.error
          = El::mpi::Recv<El::BigFloat>(from, El::mpi::COMM_WORLD);
    }
  }
}

void write_spectrum(const fs::path &output_path, const size_t &num_blocks,
                    const std::vector<Zeros> &zeros_blocks,
                    const std::vector<size_t> &block_indices)
{
  if(El::mpi::Size() == 1)
    {
      write_file(output_path, zeros_blocks);
      return;
    }

  // Synchronize zeros
  const int rank = El::mpi::Rank();

  ASSERT_EQUAL(block_indices.size(), zeros_blocks.size());

  std::vector<int> block_ranks(num_blocks, -1);
  std::map<size_t, size_t> block_index_global_to_local;
  std::vector<Zeros> zeros_all_blocks(num_blocks);

  for(size_t local_index = 0; local_index < block_indices.size();
      ++local_index)
    {
      const auto block_index = block_indices.at(local_index);
      block_ranks.at(block_index) = rank;
      block_index_global_to_local.emplace(block_index, local_index);
      zeros_all_blocks.at(block_index) = zeros_blocks.at(local_index);
    }
  El::mpi::AllReduce(block_ranks.data(), block_ranks.size(), El::mpi::MAX,
                     El::mpi::COMM_WORLD);

  for(size_t block_index = 0; block_index < num_blocks; ++block_index)
    {
      const int block_rank = block_ranks.at(block_index);
      ASSERT(block_rank >= 0, DEBUG_STRING(block_index));
      ASSERT(block_rank < El::mpi::Size(), DEBUG_STRING(block_index));

      auto &zeros_block = zeros_all_blocks.at(block_index);
      synchronize_zeros_block(zeros_block, block_rank);
    }
  write_file(output_path, zeros_all_blocks);
}
