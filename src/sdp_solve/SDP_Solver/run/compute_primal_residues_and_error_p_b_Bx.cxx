#include "sdp_solve/SDP_Solver.hxx"

// Compute the residue
//
// p[n] = dualObjective[n] - (FreeVarMatrix^T x)_n  for 0 <= n < N
//
// and the corresponding primal error max(|p_i|)

void compute_primal_residues_and_error_p_b_Bx(const Block_Info &block_info,
                                              const SDP &sdp,
                                              const Block_Vector &x,
                                              Block_Vector &primal_residue_p,
                                              El::BigFloat &primal_error)
{
  auto free_var_matrix_block(sdp.free_var_matrix.blocks.begin());
  auto x_block(x.blocks.begin());
  auto primal_residue_p_block(primal_residue_p.blocks.begin());

  El::Matrix<El::BigFloat> primal_residue_local;
  Zeros(primal_residue_local, sdp.dual_objective_b.Height(),
        sdp.dual_objective_b.Width());

  for(auto &block_index : block_info.block_indices)
    {
      El::Gemv(El::OrientationNS::TRANSPOSE, El::BigFloat(-1),
               *free_var_matrix_block, *x_block, El::BigFloat(0),
               *primal_residue_p_block);

      // The total primal error is the sum of all of the different
      // blocks.  So to prevent double counting, only add
      // dual_objective_b to one of the residue blocks.
      if(block_index == 0)
        {
          El::Axpy(El::BigFloat(1), sdp.dual_objective_b,
                   *primal_residue_p_block);
        }

      // Locally sum contributions to the primal errror
      for(int64_t row = 0; row < sdp.dual_objective_b.LocalHeight(); ++row)
        {
          int64_t global_row(sdp.dual_objective_b.GlobalRow(row));
          for(int64_t column = 0; column < sdp.dual_objective_b.LocalWidth();
              ++column)
            {
              int64_t global_column(sdp.dual_objective_b.GlobalCol(column));
              primal_residue_local(global_row, global_column)
                += primal_residue_p_block->GetLocal(row, column);
            }
        }

      ++free_var_matrix_block;
      ++x_block;
      ++primal_residue_p_block;
    }

  // Send out updates for the primal residue
  El::DistMatrix<El::BigFloat> primal_residue_dist;
  Zeros(primal_residue_dist, primal_residue_local.Height(),
        primal_residue_local.Width());

  El::BigFloat zero(0);
  for(int64_t row = 0; row < primal_residue_local.Height(); ++row)
    for(int64_t column = 0; column < primal_residue_local.Width(); ++column)
      {
        if(primal_residue_local(row, column) != zero)
          {
            primal_residue_dist.QueueUpdate(row, column,
                                            primal_residue_local(row, column));
          }
      }
  primal_residue_dist.ProcessQueues();

  // Get the max error.
  El::BigFloat local_primal_error(0);
  for(int64_t row = 0; row < primal_residue_dist.LocalHeight(); ++row)
    for(int64_t column = 0; column < primal_residue_dist.LocalWidth();
        ++column)
      {
        local_primal_error
          = std::max(local_primal_error,
                     El::Abs(primal_residue_dist.GetLocal(row, column)));
      }

  primal_error = El::mpi::AllReduce(local_primal_error, El::mpi::MAX,
                                    El::mpi::COMM_WORLD);
}
