#include "../Block_Info.hxx"

void Block_Info::compute_block_grid_mapping()
{
  El::mpi::CommGroup(El::mpi::COMM_WORLD, group);

  // Create an mpi::Group for each set of processors.  For now, one
  // group per processor.

  int rank(El::mpi::Rank(El::mpi::COMM_WORLD));
  El::mpi::Incl(group, 1, &rank, group);
  El::mpi::Create(El::mpi::COMM_WORLD, group, comm);

  size_t num_procs(El::mpi::Size(El::mpi::COMM_WORLD));

  // FIXME: Super simple round-robin assigning of MPI processes per
  // block.  There is no attempt at sophisticated balancing.

  const size_t num_blocks(dimensions.size());
  for(size_t jj = 0; jj < num_blocks; ++jj)
    {
      int processor(jj % num_procs);
      if(processor == El::mpi::Rank())
        {
          block_indices.push_back(jj);
        }
    }
}
