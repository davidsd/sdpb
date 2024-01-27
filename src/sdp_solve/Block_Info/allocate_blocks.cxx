#include "../Block_Info.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/block_mapping/compute_block_grid_mapping.hxx"
#include "sdpb_util/block_mapping/create_mpi_block_mapping_groups.hxx"

void Block_Info::allocate_blocks(const Environment &env,
                                 const std::vector<Block_Cost> &block_costs,
                                 const size_t &proc_granularity,
                                 const Verbosity &verbosity)
{
  const size_t num_procs = El::mpi::Size(El::mpi::COMM_WORLD);
  const size_t num_nodes = env.num_nodes();
  const size_t procs_per_node = env.comm_shared_mem.Size();
  ASSERT(procs_per_node * num_nodes == num_procs,
         "Incompatible number of MPI processes and processes per node "
         "for node=",
         env.node_index(),
         ". Each node should have the same number of processes."
         "\n\tMPI processes: ",
         num_procs, "\n\tprocsPerNode: ", procs_per_node,
         "\n\tnum_nodes: ", num_nodes);
  ASSERT(procs_per_node % proc_granularity == 0,
         "Incompatible number of processes per node and process granularity.  "
         "procGranularity mush evenly divide procsPerNode:\n\tprocsPerNode: ",
         procs_per_node, "\n\tprocGranularity: ", proc_granularity);

  std::vector<std::vector<Block_Map>> mapping(compute_block_grid_mapping(
    procs_per_node / proc_granularity, num_nodes, block_costs));

  for(auto &block_vector : mapping)
    for(auto &block_map : block_vector)
      {
        block_map.num_procs *= proc_granularity;
      }

  if(verbosity >= Verbosity::debug && El::mpi::Rank() == 0)
    {
      std::stringstream ss;
      ss << "Block Grid Mapping\n"
         << "Node\tNum Procs\tCost (Per Proc)\t\tBlock List\n"
         << "==========================================================\n";
      for(size_t node = 0; node < mapping.size(); ++node)
        {
          for(auto &m : mapping[node])
            {
              ss << node << "\t" << m.num_procs << "\t\t"
                 << m.cost / static_cast<double>(m.num_procs) << "\t\t\t{";
              for(size_t ii = 0; ii < m.block_indices.size(); ++ii)
                {
                  if(ii != 0)
                    {
                      ss << ", ";
                    }
                  ss << m.block_indices[ii] << "("
                     << dimensions[m.block_indices[ii]] << ","
                     << num_points[m.block_indices[ii]] << ")";
                }
              ss << "}\n";
            }
          ss << "\n";
        }
      El::Output(ss.str());
    }

  const auto &node_comm = env.comm_shared_mem;
  const int node_index = env.node_index();

  create_mpi_block_mapping_groups(mapping, node_comm, node_index,
                                  mpi_group.value, mpi_comm.value,
                                  block_indices);

  ASSERT(!block_indices.empty(),
         "No SDP blocks were assigned to rank=", El::mpi::Rank(),
         ". node=", node_index, " node_rank=", node_comm.Rank());
}
