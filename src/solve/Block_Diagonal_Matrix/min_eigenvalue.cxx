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
      Real min_block_lambda
        = min_eigenvalue(A.blocks[b], workspace[b], eigenvalues[b]);
      // ensure only one thread modifies lambdaMin at a time
      {
        lambda_min = min(lambda_min, min_block_lambda);
      }
    }

  return lambda_min;
}
