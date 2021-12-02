#pragma once

#include "sdp_solve/Block_Info.hxx"
#include "sdp_solve/Write_Solution.hxx"
#include "sdp_solve/read_text_block.hxx"
#include "dynamical_solve/Dynamical_Solver.hxx"

El::BigFloat dot(const Block_Vector &A, const Block_Vector &B);
