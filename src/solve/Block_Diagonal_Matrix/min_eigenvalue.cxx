#include "../Block_Diagonal_Matrix.hxx"

// Minimum eigenvalue of A, via the QR method
// Inputs:
// A           : symmetric Block_Diagonal_Matrix
//
El::BigFloat min_eigenvalue(Block_Diagonal_Matrix &A)
{
  El::BigFloat lambda_min(El::limits::Max<El::BigFloat>());

  El::DistMatrix<El::BigFloat> eigenvalues;
  for(auto &block : A.blocks_elemental)
    {
      /// FIXME: Maybe use regular double precision here?
      El::HermitianEig(El::UpperOrLowerNS::LOWER, block, eigenvalues);
      lambda_min = El::Min(lambda_min, El::Min(eigenvalues));
    }
  return lambda_min;
}
