#include "../Block_Diagonal_Matrix.hxx"

// Minimum eigenvalue of A.  A is assumed to be symmetric.

// Annoyingly, El::HermitianEig modifies 'block'.  It is OK, because
// it is only called from step_length(), which passes in a temporary.
// Still ugly.
El::BigFloat min_eigenvalue(Block_Diagonal_Matrix &A)
{
  El::BigFloat local_min(El::limits::Max<El::BigFloat>());

  for(auto &block : A.blocks)
    {
      El::DistMatrix<El::BigFloat, El::VR, El::STAR> eigenvalues(block.Grid());
      El::HermitianEig(El::UpperOrLowerNS::LOWER, block, eigenvalues);
      local_min = El::Min(local_min, El::Min(eigenvalues));
    }
  return El::mpi::AllReduce(local_min, El::mpi::MIN, El::mpi::COMM_WORLD);
}
