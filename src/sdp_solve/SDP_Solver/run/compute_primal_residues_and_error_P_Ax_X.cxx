#include "constraint_matrix_weighted_sum.hxx"

// PrimalResidues = \sum_p A_p x[p] - X

void compute_primal_residues_and_error_P_Ax_X(
  const Block_Info &block_info, const SDP &sdp, const Block_Vector &x,
  const Block_Diagonal_Matrix &X, Block_Diagonal_Matrix &primal_residues,
  El::BigFloat &primal_error, Timers &timers)
{
  Scoped_Timer primal_residues_timer(timers, "computePrimalResidues");
  constraint_matrix_weighted_sum(block_info, sdp, x, primal_residues);
  primal_residues -= X;
  primal_error = primal_residues.max_abs();
}
