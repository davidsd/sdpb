#include "../Block_Diagonal_Matrix.hxx"

// Minimum eigenvalue of A, via the QR method
// Inputs:
// A           : symmetric Block_Diagonal_Matrix
// eigenvalues : vector<Vector> of length A.blocks.size()
// workspace   : vector<Vector> of length A.blocks.size()
//
El::BigFloat min_eigenvalue(Block_Diagonal_Matrix &A)
{
  El::BigFloat lambda_min(El::limits::Max<El::BigFloat>());
  for(auto &block : A.blocks_elemental)
    {
      /// FIXME: Maybe use regular double precision here?
      El::DistMatrix<El::BigFloat> eigenvalues;
      El::HermitianEig(El::UpperOrLowerNS::LOWER, block, eigenvalues);
      lambda_min = El::Min(lambda_min, El::Min(eigenvalues));
    }
  return lambda_min;
}
