#include "constraint_matrix_weighted_sum.hxx"
#include "sdpa_solve/SDP.hxx"

namespace Sdpb::Sdpa
{
  // PrimalResidues = \sum_p F_p x[p] - F_0 - X
  void
  compute_primal_residues_and_error(const SDP &sdp,
                                    const Primal_Dist_Vector &x,
                                    const Block_Diagonal_Matrix &X,
                                    Block_Diagonal_Matrix &primal_residues,
                                    El::BigFloat &primal_error, Timers &timers)
  {
    Scoped_Timer primal_residues_timer(timers, "computePrimalResidues");
    constraint_matrix_weighted_sum(sdp, x, primal_residues);
    primal_residues -= sdp.sdp_block_F_0;
    primal_residues -= X;
    primal_error = primal_residues.max_abs();
  }
}
