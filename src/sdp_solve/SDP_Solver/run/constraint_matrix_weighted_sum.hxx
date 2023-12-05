#pragma once

#include "sdp_solve/SDP_Solver.hxx"

void constraint_matrix_weighted_sum(const Block_Info &block_info,
                                    const SDP &sdp, const Block_Vector &a,
                                    Block_Diagonal_Matrix &Result);
