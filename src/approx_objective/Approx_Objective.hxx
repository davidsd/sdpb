#pragma once

#include "sdp_solve/sdp_solve.hxx"

struct Approx_Objective
{
  El::BigFloat objective, d_objective, dd_objective;
  Approx_Objective(const Block_Info &block_info, const SDP &sdp,
                   const SDP &d_sdp, const Block_Vector &x,
                   const Block_Vector &y,
                   const Block_Diagonal_Matrix &schur_complement_cholesky,
                   const Block_Matrix &schur_off_diagonal,
                   const El::DistMatrix<El::BigFloat> &Q);

  Approx_Objective(const SDP &sdp, const SDP &d_sdp, const Block_Vector &x,
                   const Block_Vector &y);
};
