#include "BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdp_solve/Block_Matrix.hxx"
#include "sdp_solve/Block_Info.hxx"
#include "sdp_solve/SDP.hxx"

inline BigInt_Shared_Memory_Syrk_Context
initialize_bigint_syrk_context(const Environment &env,
                               const Block_Info &block_info, const SDP &sdp,
                               bool debug)
{
  int block_width = sdp.dual_objective_b.Height(); // = N

  // Communicator for all ranks on a node
  const auto &shared_memory_comm = env.comm_shared_mem;

  auto num_ranks_per_node = El::mpi::Size(shared_memory_comm);

  // Collect block heights and block indices from all ranks in shared_memory_comm
  std::map<int, std::vector<size_t>> rank_to_global_block_indices;
  std::map<int, std::vector<int>> rank_to_block_heights;
  for(int rank = 0; rank < num_ranks_per_node; ++rank)
    {
      size_t num_blocks_in_rank = sdp.free_var_matrix.blocks.size();
      assert(num_blocks_in_rank == block_info.block_indices.size());
      El::mpi::Broadcast(num_blocks_in_rank, rank, shared_memory_comm);

      rank_to_global_block_indices.emplace(
        rank, std::vector<size_t>(num_blocks_in_rank));
      rank_to_block_heights.emplace(rank,
                                    std::vector<int>(num_blocks_in_rank));

      if(num_blocks_in_rank == 0)
        continue;
      if(El::mpi::Rank(shared_memory_comm) == rank)
        {
          rank_to_global_block_indices[rank] = block_info.block_indices;
          for(size_t block_index = 0; block_index < num_blocks_in_rank;
              ++block_index)
            {
              rank_to_block_heights[rank][block_index]
                = sdp.free_var_matrix.blocks[block_index].Height();
            }
        }
      assert(num_blocks_in_rank == rank_to_global_block_indices[rank].size());
      assert(num_blocks_in_rank == rank_to_block_heights[rank].size());

      El::mpi::Broadcast(rank_to_global_block_indices[rank].data(),
                         num_blocks_in_rank, rank, shared_memory_comm);
      El::mpi::Broadcast(rank_to_block_heights[rank].data(),
                         num_blocks_in_rank, rank, shared_memory_comm);
    }

  std::map<size_t, int> global_block_index_to_height;
  std::set<size_t> global_block_indices;
  for(int rank = 0; rank < num_ranks_per_node; ++rank)
    {
      for(size_t i = 0; i < rank_to_global_block_indices[rank].size(); ++i)
        {
          auto global_index = rank_to_global_block_indices[rank][i];
          global_block_indices.insert(global_index);

          auto height = rank_to_block_heights[rank][i];
          global_block_index_to_height.emplace(global_index, height);
        }
    }

  std::map<size_t, size_t> block_index_global_to_shmem;
  std::vector<int> shmem_block_index_to_height(global_block_indices.size());
  size_t curr_shmem_block_index = 0;
  for(size_t global_index : global_block_indices)
    {
      block_index_global_to_shmem.emplace(global_index,
                                          curr_shmem_block_index);
      shmem_block_index_to_height[curr_shmem_block_index]
        = global_block_index_to_height[global_index];
      curr_shmem_block_index++;
    }

  std::vector<size_t> block_index_local_to_shmem(
    block_info.block_indices.size());
  for(size_t local_index = 0; local_index < block_index_local_to_shmem.size();
      ++local_index)
    {
      size_t global_index = block_info.block_indices[local_index];
      size_t shmem_index = block_index_global_to_shmem[global_index];
      block_index_local_to_shmem[local_index] = shmem_index;
    }

  // Print block indices and sizes
  if(debug && El::mpi::Rank(shared_memory_comm) == 0)
    {
      std::ostringstream os;
      El::BuildStream(
        os, "initialize_bigint_syrk_context, node=", env.node_index(), "\n");

      // global block indices on the node
      std::vector<int> block_index_shmem_to_global(
        block_index_global_to_shmem.size());
      for(const auto &[global_index, shmem_index] :
          block_index_global_to_shmem)
        {
          block_index_shmem_to_global.at(shmem_index) = global_index;
        }
      El::BuildStream(os, "Number of blocks on the node: ",
                      block_index_shmem_to_global.size(), "\n");
      El::Print(block_index_shmem_to_global,
                "Block indices on the node: ", ", ", os);
      os << "\n";

      El::Print(shmem_block_index_to_height, "Blocks heights: ", ", ", os);
      os << "\n";

      auto total_height
        = std::accumulate(shmem_block_index_to_height.begin(),
                          shmem_block_index_to_height.end(), 0);
      El::BuildStream(os, "Total block height: ", total_height, "\n");
      El::BuildStream(os, "Block width: ", block_width, "\n");
      El::BuildStream(os, "Total block elements: ", total_height * block_width,
                      "\n");
      El::BuildStream(os, "Total Q elements: ", block_width * block_width,
                      "\n");

      El::Output(os.str());
    }

  return BigInt_Shared_Memory_Syrk_Context(
    shared_memory_comm, El::gmp::Precision(), shmem_block_index_to_height,
    block_width, block_index_local_to_shmem, debug);
}
