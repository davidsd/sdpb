#include "../Block_Diagonal_Matrix.hxx"

// Minimum eigenvalue of A, via the QR method
// Inputs:
// A           : symmetric Block_Diagonal_Matrix
// eigenvalues : vector<Vector> of length A.blocks.size()
// workspace   : vector<Vector> of length A.blocks.size()
//
Real min_eigenvalue(Block_Diagonal_Matrix &A, std::vector<Vector> &workspace,
                    std::vector<Vector> &eigenvalues)
{
  assert(A.blocks.size() == eigenvalues.size());
  assert(A.blocks.size() == workspace.size());

  // TODO(davidsd): get rid of this hack
  Real lambda_min = 1e100; // we really want lambdaMin = infinity here
  for(size_t b = 0; b < A.blocks.size(); b++)
    {
      lambda_min = min(
        lambda_min, min_eigenvalue(A.blocks[b], workspace[b], eigenvalues[b]));
    }

  return lambda_min;
}

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
