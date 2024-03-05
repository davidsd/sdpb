#include "sdp_solve/Block_Diagonal_Matrix.hxx"

// X := ACholesky^{-T} ACholesky^{-1} X = A^{-1} X
void cholesky_solve(const Block_Diagonal_Matrix &ACholesky,
                    Block_Diagonal_Matrix &X)
{
  for(size_t b = 0; b < X.blocks.size(); b++)
    {
      El::cholesky::SolveAfter(El::UpperOrLowerNS::LOWER,
                               El::OrientationNS::NORMAL, ACholesky.blocks[b],
                               X.blocks[b]);
    }
}
