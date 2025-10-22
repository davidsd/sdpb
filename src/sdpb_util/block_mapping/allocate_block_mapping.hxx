#pragma once

#include "sdpb_util/assert.hxx"
#include "sdpb_util/Environment.hxx"
#include "sdpb_util/Verbosity.hxx"
#include "compute_block_grid_mapping.hxx"
#include "create_mpi_block_mapping_groups.hxx"

inline std::vector<std::vector<Block_Map>> allocate_block_mapping(
  const Environment &env, const std::vector<Block_Cost> &block_costs,
  const size_t &proc_granularity,
  const std::function<std::ostream &(std::ostream &, size_t)> &print_block,
  const Verbosity &verbosity, El::mpi::Group &mpi_group,
  El::mpi::Comm &mpi_group_comm, std::vector<size_t> &block_indices)
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
                  print_block(ss, m.block_indices[ii]);
                }
              ss << "}\n";
            }
          ss << "\n";
        }
      El::Output(ss.str());
    }

  if(El::mpi::Rank() == 0)
    {
      for(size_t node = 0; node < mapping.size(); ++node)
        {
          const auto &node_mapping = mapping[node];
          for(size_t group = 0; group < node_mapping.size(); ++group)
            {
              const auto &group_mapping = node_mapping[group];
              ASSERT(!group_mapping.block_indices.empty(),
                     "No SDP blocks were assigned to node=", node,
                     ", group=", group);
            }
        }
    }

  const auto &node_comm = env.comm_shared_mem;
  const int node_index = env.node_index();

  create_mpi_block_mapping_groups(mapping, node_comm, node_index, mpi_group,
                                  mpi_group_comm, block_indices);

  ASSERT(!block_indices.empty(),
         "No SDP blocks were assigned to rank=", El::mpi::Rank(),
         ". node=", node_index, " node_rank=", node_comm.Rank());
  return mapping;
}
