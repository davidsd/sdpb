#pragma once

#include "Block_Diagonal_Matrix.hxx"
#include "Block_Vector.hxx"

// v := L^{-T} v, where L is lower-triangular
void lower_triangular_transpose_solve(const Block_Diagonal_Matrix &L, Block_Vector &v);
