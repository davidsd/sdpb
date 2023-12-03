#include "sdp_solve/sdp_solve.hxx"

void Axpy(const El::BigFloat &alpha, const SDP &new_sdp, SDP &delta_sdp)
{
  // Assume that the new and old bilinear bases are the same
  for(size_t block(0); block!=delta_sdp.free_var_matrix.blocks.size(); ++block)
  {
    El::Axpy(alpha, new_sdp.free_var_matrix.blocks.at(block),
             delta_sdp.free_var_matrix.blocks.at(block));
    El::Axpy(alpha, new_sdp.primal_objective_c.blocks.at(block),
             delta_sdp.primal_objective_c.blocks.at(block));
  }
  El::Axpy(alpha, new_sdp.dual_objective_b, delta_sdp.dual_objective_b);
  delta_sdp.objective_const += alpha*new_sdp.objective_const;
}
