#include "../Block_Diagonal_Matrix.hxx"

// Minimum eigenvalue of A.  A is assumed to be symmetric.

El::BigFloat min_eigenvalue(Block_Diagonal_Matrix &A)
{
  El::BigFloat lambda_min(El::limits::Max<El::BigFloat>());

  for(auto block=A.blocks_elemental.begin(); block!=A.blocks_elemental.end(); ++block)
    {
      El::DistMatrix<El::BigFloat,El::VR,El::STAR> eigenvalues( block->Grid() );

      El::DistMatrix<El::BigFloat> block_copy(block->Height(),block->Width());
      block_copy=*block;

      El::HermitianEig(El::UpperOrLowerNS::LOWER, block_copy, eigenvalues);
      lambda_min = El::Min(lambda_min, El::Min(eigenvalues));
    }
  return lambda_min;
}
