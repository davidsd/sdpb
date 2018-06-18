#include "Block_Diagonal_Matrix.hxx"
#include "Block_Vector.hxx"

// v := L^{-T} v, where L is lower-triangular
void lower_triangular_transpose_solve(const Block_Diagonal_Matrix &L,
                                      Block_Vector &v)
{
  for(size_t b = 0; b < L.blocks.size(); b++)
    {
      El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
               El::OrientationNS::TRANSPOSE, El::UnitOrNonUnit::NON_UNIT,
               El::BigFloat(1), L.blocks[b], v.blocks[b]);
    }
}
