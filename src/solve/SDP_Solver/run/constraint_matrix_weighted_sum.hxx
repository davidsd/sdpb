#pragma once

#include "../../SDP_Solver.hxx"

void constraint_matrix_weighted_sum(const SDP &sdp, const Vector &a,
                                    Block_Diagonal_Matrix &Result);
