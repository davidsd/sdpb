#include "../../../SDP_Solver.hxx"

// Compute the residue
//
// p[n] = dualObjective[n] - (FreeVarMatrix^T x)_n  for 0 <= n < N
//
// and the corresponding primal error max(|p_i|)

void compute_primal_residues_and_error_p(const Block_Info &block_info,
                                         const SDP &sdp, const Block_Vector &x,
                                         Block_Vector &primal_residue_p,
                                         El::BigFloat &primal_error)
{
  auto free_var_matrix_block(sdp.free_var_matrix.blocks.begin());
  auto x_block(x.blocks.begin());
  auto primal_residue_p_block(primal_residue_p.blocks.begin());

  El::BigFloat local_primal_error(primal_error);
  for(auto &block_index : block_info.block_indices)
    {
      El::Gemv(El::OrientationNS::TRANSPOSE, El::BigFloat(-1),
               *free_var_matrix_block, *x_block, El::BigFloat(0),
               *primal_residue_p_block);

      // Contribution to the primal errror |p_i|
      for(int64_t row = 0; row < sdp.dual_objective_b.LocalHeight(); ++row)
        {
          for(int64_t column = 0; column < sdp.dual_objective_b.LocalWidth();
              ++column)
            {
              local_primal_error = std::max(
                local_primal_error,
                El::Abs(sdp.dual_objective_b.GetLocal(row, column)
                        + primal_residue_p_block->GetLocal(row, column)));
            }
        }
      if(block_index == 0)
        {
          El::Axpy(El::BigFloat(1), sdp.dual_objective_b,
                   *primal_residue_p_block);
        }
      ++free_var_matrix_block;
      ++x_block;
      ++primal_residue_p_block;
    }
  primal_error = El::mpi::AllReduce(local_primal_error, El::mpi::MAX,
                                    El::mpi::COMM_WORLD);
}
