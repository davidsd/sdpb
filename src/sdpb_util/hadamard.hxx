#pragma once

#include "assert.hxx"

#include <El.hpp>

// Elementwise multiplication of DistMatrix<STAR,STAR> (i.e. each rank holds a copy of a whole matrix)
// and DistMatrix with arbitrary distribution.
template <class Ring, El::Dist U, El::Dist V>
void hadamard(const El::DistMatrix<Ring, El::STAR, El::STAR> &A,
              const El::DistMatrix<Ring, U, V> &B,
              El::DistMatrix<Ring, U, V> &C)
{
  ASSERT(A.Grid() == B.Grid());
  ASSERT(A.Grid() == C.Grid());
  ASSERT_EQUAL(A.Height(), B.Height());
  ASSERT_EQUAL(A.Width(), B.Width());

  C.Resize(A.Height(), A.Width());

  for(int iLoc = 0; iLoc < C.LocalHeight(); ++iLoc)
    for(int jLoc = 0; jLoc < C.LocalWidth(); ++jLoc)
      {
        const int i = C.GlobalRow(iLoc);
        const int j = B.GlobalCol(jLoc);

        const auto &a = A.LockedMatrix()(i, j);
        const auto &b = B.LockedMatrix()(iLoc, jLoc);
        C.SetLocal(iLoc, jLoc, a * b);
      }
}

// This combination is already implemented in Elemental
template <class Ring>
void hadamard(const El::DistMatrix<Ring, El::STAR, El::STAR> &A,
              const El::DistMatrix<Ring, El::STAR, El::STAR> &B,
              El::DistMatrix<Ring, El::STAR, El::STAR> &C)
{
  El::Hadamard(A, B, C);
}
