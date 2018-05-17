#include "../Block_Diagonal_Matrix.hxx"

// A := L^{-1} A L^{-T}
void lower_triangular_inverse_congruence(Block_Diagonal_Matrix &A,
                                         Block_Diagonal_Matrix &L)
{
  for(size_t b = 0; b < A.blocks.size(); b++)
    {
      lower_triangular_inverse_congruence(A.blocks[b], L.blocks[b]);

      El::Trsm(El::LeftOrRight::RIGHT, El::UpperOrLowerNS::LOWER,
               El::Orientation::TRANSPOSE, El::UnitOrNonUnit::NON_UNIT,
               El::BigFloat(1), L.blocks_elemental[b], A.blocks_elemental[b]);
      El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
               El::Orientation::NORMAL, El::UnitOrNonUnit::NON_UNIT,
               El::BigFloat(1), L.blocks_elemental[b], A.blocks_elemental[b]);
    }
}
