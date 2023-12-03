#pragma once

#include "sdp_solve/Block_Diagonal_Matrix.hxx"

// B := L^{-1} B, where L is the result of a previous cholesky
// factorization.  Note that this is different from computing the solution to
// A B == (L L^T) B

template <class T>
void lower_triangular_solve(const Block_Diagonal_Matrix &L_cholesky, T &B)
{
  for(size_t block = 0; block < L_cholesky.blocks.size(); block++)
    {
      El::Trsm(El::LeftOrRightNS::LEFT, El::UpperOrLowerNS::LOWER,
               El::OrientationNS::NORMAL, El::UnitOrNonUnitNS::NON_UNIT,
               El::BigFloat(1), L_cholesky.blocks[block], B.blocks[block]);
    }
}
