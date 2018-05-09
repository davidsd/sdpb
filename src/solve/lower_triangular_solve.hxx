#pragma once

#include "Block_Diagonal_Matrix.hxx"
#include "Block_Matrix.hxx"

// B := L^{-1} B, where L is the result of a previous cholesky
// factorization.  Note that this is different from computing the solution to
// A B == (L L^T) B
void lower_triangular_solve(const Block_Diagonal_Matrix &L, Block_Matrix &B);
