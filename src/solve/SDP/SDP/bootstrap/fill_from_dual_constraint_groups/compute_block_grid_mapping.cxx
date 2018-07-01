#include "../../../../SDP.hxx"

void compute_block_grid_mapping(const size_t &num_blocks,
                                std::vector<size_t> &block_indices)
{
  El::mpi::Group default_group;
  El::mpi::CommGroup(El::mpi::COMM_WORLD, default_group);

  // Create an mpi::Group for each set of processors.  For now, one
  // group per processor.
  //
  // El::mpi::Group group(default_group);
  // El::mpi::Incl(default_group, 1, &processor, group);
  // El::mpi::Comm comm(El::mpi::COMM_WORLD);
  // El::mpi::Create(El::mpi::COMM_WORLD, group, comm);

  size_t num_procs(El::mpi::Size(default_group));

  // FIXME: Super simple round-robin assigning of MPI processes per
  // block.  There is no attempt at sophisticated balancing.

  for(size_t jj = 0; jj < num_blocks; ++jj)
    {
      int processor(jj % num_procs);
      if(processor == El::mpi::Rank())
        {
          block_indices.push_back(jj);
        }
    }
}
