#include "sdp_solve/Block_Diagonal_Matrix.hxx"

// A := L^{-1} A L^{-T}
void lower_triangular_inverse_congruence(const Block_Diagonal_Matrix &L,
                                         Block_Diagonal_Matrix &A)
{
  for(size_t b = 0; b < A.blocks.size(); b++)
    {
      El::Trsm(El::LeftOrRight::RIGHT, El::UpperOrLowerNS::LOWER,
               El::Orientation::TRANSPOSE, El::UnitOrNonUnit::NON_UNIT,
               El::BigFloat(1), L.blocks[b], A.blocks[b]);
      El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
               El::Orientation::NORMAL, El::UnitOrNonUnit::NON_UNIT,
               El::BigFloat(1), L.blocks[b], A.blocks[b]);
    }
}
