#include "Block_Diagonal_Matrix.hxx"
#include "Block_Matrix.hxx"

// B := L^{-1} B, where L is the result of a previous cholesky
// factorization.  Note that this is different from computing the solution to
// A B == (L L^T) B
void lower_triangular_solve(const Block_Diagonal_Matrix &L, Block_Matrix &B)
{
  for(size_t b = 0; b < L.blocks.size(); b++)
    {
      El::Trsm(El::LeftOrRightNS::LEFT, El::UpperOrLowerNS::LOWER,
               El::OrientationNS::NORMAL, El::UnitOrNonUnitNS::NON_UNIT,
               El::BigFloat(1), L.blocks_elemental[b], B.blocks[b]);
    }
}
