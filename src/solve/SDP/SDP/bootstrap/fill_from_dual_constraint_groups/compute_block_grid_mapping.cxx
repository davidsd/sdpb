#include "../../../../SDP.hxx"

void compute_block_grid_mapping(const size_t &num_blocks,
                                std::list<El::Grid> &mapping)
{
  El::mpi::Group default_group;
  El::mpi::CommGroup(El::mpi::COMM_WORLD, default_group);

  // FIXME: Super simple round-robin assigning of MPI processes per
  // block.  There is no attempt at sophisticated balancing.

  size_t num_procs(El::mpi::Size(default_group));

  for(size_t jj = 0; jj < num_blocks; ++jj)
    {
      // El::mpi::Group sub_group;
      // int processor(jj % num_procs);
      // El::mpi::Incl(default_group, 1, &processor, sub_group);
      // mapping.emplace_back(El::mpi::COMM_WORLD, sub_group, 1);
      mapping.emplace_back(El::mpi::COMM_WORLD);
    }
}
