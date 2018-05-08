#pragma once

#include "../../SDP_Solver.hxx"

void constraint_matrix_weighted_sum(const SDP &sdp, const Vector &a,
                                    Block_Diagonal_Matrix &Result);

void constraint_matrix_weighted_sum(const SDP &sdp, const Block_Matrix &a,
                                    Block_Diagonal_Matrix &Result);
