#include "Block_Diagonal_Matrix.hxx"
#include "Block_Vector.hxx"

// v := L^{-T} v, where L is lower-triangular
void lower_triangular_transpose_solve(const Block_Diagonal_Matrix &L,
                                      Block_Vector &v)
{
  for(size_t b = 0; b < L.blocks_elemental.size(); b++)
    {
      El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
               El::OrientationNS::TRANSPOSE, El::UnitOrNonUnit::NON_UNIT,
               El::BigFloat(1), L.blocks_elemental[b], v.blocks[b]);
    }
}
