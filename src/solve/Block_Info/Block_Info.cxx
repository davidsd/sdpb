#include "../Block_Info.hxx"

#include "../../compute_block_grid_mapping.hxx"

Block_Info::Block_Info(const boost::filesystem::path &sdp_directory,
                       const size_t &procs_per_node)
{
  boost::filesystem::ifstream block_stream(sdp_directory / "blocks");
  if(!block_stream.good())
    {
      throw std::runtime_error("Could not open '"
                               + (sdp_directory / "blocks").string() + "'");
    }
  read_vector(block_stream, dimensions);
  read_vector(block_stream, degrees);
  read_vector(block_stream, schur_block_sizes);
  read_vector(block_stream, psd_matrix_block_sizes);
  read_vector(block_stream, bilinear_pairing_block_sizes);

  const size_t num_procs(El::mpi::Size(El::mpi::COMM_WORLD));
  std::vector<Block_Cost> block_costs;
  // Assume that the cost of a block is about N^3
  for(size_t ii = 0; ii < schur_block_sizes.size(); ++ii)
    {
      block_costs.emplace_back(schur_block_sizes[ii] * schur_block_sizes[ii]
                                 * schur_block_sizes[ii],
                               ii);
    }
  // Reverse sort, with largest first
  std::sort(block_costs.rbegin(), block_costs.rend());
  if(num_procs % procs_per_node != 0)
    {
      throw std::runtime_error(
        "Incompatible number of processes and processes per node.  "
        "procs_per_node must evenly divide num_procs:\n\tnum_procs: "
        + std::to_string(num_procs)
        + "\n\tprocs_per_node: " + std::to_string(procs_per_node));
    }
  const size_t num_nodes(num_procs / procs_per_node);
  std::vector<std::vector<Block_Map>> mapping(
    compute_block_grid_mapping(procs_per_node, num_nodes, block_costs));

  // Create an mpi::Group for each set of processors.

  El::mpi::Group default_mpi_group;
  El::mpi::CommGroup(El::mpi::COMM_WORLD, default_mpi_group);

  int rank(El::mpi::Rank(El::mpi::COMM_WORLD));
  int current_rank(0);
  for(auto &block_vector : mapping)
    {
      for(auto &block_map : block_vector)
        {
          current_rank += block_map.num_procs;
          if(current_rank > rank)
            {
              block_indices = block_map.block_indices;
              break;
            }
        }
      if(!block_indices.empty())
        {
          break;
        }
    }
  {
    // MPI wants 'int' arrays, not 'size_t' arrays
    std::vector<int> int_indices;
    std::copy(block_indices.begin(), block_indices.end(), int_indices.begin());
    El::mpi::Incl(default_mpi_group, int_indices.size(), int_indices.data(),
                  mpi_group);
  }
  El::mpi::Create(El::mpi::COMM_WORLD, mpi_group, mpi_comm);
}
