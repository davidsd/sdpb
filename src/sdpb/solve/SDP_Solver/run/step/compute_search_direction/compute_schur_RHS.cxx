#include "../../../../SDP_Solver.hxx"

// Compute the vectors dx and dy on the right-hand side of the Schur
// complement equation:
//
// {{S, -B}, {B^T, 0}} . {dx, dy} = {dx, dy}
//
// where S = SchurComplement and B = FreeVarMatrix.  Specifically,
//
// dx[p] = -dual_residues[p] - Tr(A_p Z)              for 0 <= p < P
// dy[n] = dualObjective[n] - (FreeVarMatrix^T x)_n  for 0 <= n < N
//
// where P = primalObjective.size(), N = dualObjective.size()
//
// Inputs:
// - sdp, an SDP
// - dual_residues, a Vector of length P
// - Z = X^{-1} (PrimalResidues Y - R)
// - x, a vector of length P
// Outputs:
// - dx, a Vector of length P
// - dy, a Vector of length N
//

void compute_schur_RHS(const Block_Info &block_info, const SDP &sdp,
                       const Block_Vector &dual_residues,
                       const Block_Diagonal_Matrix &Z, const Block_Vector &x,
                       Block_Vector &dx, Block_Vector &dy)
{
  auto free_var_matrix_block(sdp.free_var_matrix.blocks.begin());
  auto dual_residues_block(dual_residues.blocks.begin());
  auto x_block(x.blocks.begin());
  auto dx_block(dx.blocks.begin());
  auto dy_block(dy.blocks.begin());

  auto Z_block(Z.blocks.begin());
  auto bilinear_bases_block(sdp.bilinear_bases_dist.begin());
  for(auto &block_index : block_info.block_indices)
    {
      // dx = -dual_residues
      *dx_block = *dual_residues_block;
      *dx_block *= -1;
      const size_t dx_block_size(block_info.degrees[block_index] + 1);

      // dx[p] -= Tr(A_p Z)
      // Not sure whether it is better to first loop over blocks in
      // the result or over sub-blocks in Z
      for(size_t parity = 0; parity < 2; ++parity)
        {
          const size_t Z_block_size(bilinear_bases_block->Height());
          El::DistMatrix<El::BigFloat> ones(Z_block->Grid());
          El::Ones(ones, Z_block_size, 1);

          for(size_t column_block = 0;
              column_block < block_info.dimensions[block_index];
              ++column_block)
            for(size_t row_block = 0; row_block <= column_block; ++row_block)
              {
                size_t column_offset(column_block * Z_block_size),
                  row_offset(row_block * Z_block_size);

                El::DistMatrix<El::BigFloat> Z_sub_block(
                  El::LockedView(*Z_block, row_offset, column_offset,
                                 Z_block_size, Z_block_size)),
                  Z_times_q(Z_block_size, dx_block_size, Z_block->Grid());
                El::Zero(Z_times_q);
                El::DistMatrix<El::BigFloat> q_Z_q(Z_times_q);
                El::Zero(q_Z_q);

                El::Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL,
                         El::BigFloat(1), Z_sub_block,
                         *bilinear_bases_block, El::BigFloat(0),
                         Z_times_q);

                El::Hadamard(Z_times_q, *bilinear_bases_block, q_Z_q);

                const size_t dx_row_offset(
                  ((column_block * (column_block + 1)) / 2 + row_block)
                  * dx_block_size);
                El::DistMatrix<El::BigFloat> dx_sub_block(
                  El::View(*dx_block, dx_row_offset, 0, dx_block_size, 1));

                El::Gemv(El::Orientation::TRANSPOSE, El::BigFloat(-1), q_Z_q,
                         ones, El::BigFloat(1), dx_sub_block);
              }
          ++Z_block;
          ++bilinear_bases_block;
        }

      // dy = dualObjective - (FreeVarMatrix^T x)
      El::Gemv(El::OrientationNS::TRANSPOSE, El::BigFloat(-1),
               *free_var_matrix_block, *x_block, El::BigFloat(0), *dy_block);

      // The dy.blocks each only have part of the entire dy vector.
      // It will later all get summed together, so add
      // dual_objective_b to only the first block.
      if(block_index == 0)
        {
          El::Axpy(El::BigFloat(1), sdp.dual_objective_b, *dy_block);
        }
      ++free_var_matrix_block;
      ++dual_residues_block;
      ++x_block;
      ++dx_block;
      ++dy_block;
    }
}
