#include "../../../SDP_Solver.hxx"

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

void compute_schur_RHS(const SDP &sdp, const Block_Vector &dual_residues,
                       const Block_Diagonal_Matrix &Z, const Block_Vector &x,
                       Block_Vector &dx, Block_Vector &dy)
{
  for(size_t jj = 0; jj < sdp.dimensions.size(); ++jj)
    {
      // dx = -dual_residues
      dx.blocks[jj] = dual_residues.blocks[jj];
      dx.blocks[jj] *= -1;
      const size_t dx_block_size(sdp.degrees[jj] + 1);

      // dx[p] -= Tr(A_p Z)
      // Not sure whether it is better to first loop over blocks in
      // the result or over sub-blocks in Z
      for(size_t block_index = 2 * jj; block_index < 2 * jj + 2; ++block_index)
        {
          const size_t Z_block_size(
            sdp.bilinear_bases_dist[block_index].Height());
          El::DistMatrix<El::BigFloat> ones(Z.blocks[block_index].Grid());
          El::Ones(ones, Z_block_size, 1);

          for(size_t column_block = 0; column_block < sdp.dimensions[jj];
              ++column_block)
            for(size_t row_block = 0; row_block <= column_block; ++row_block)
              {
                size_t column_offset(column_block * Z_block_size),
                  row_offset(row_block * Z_block_size);

                El::DistMatrix<El::BigFloat> Z_sub_block(
                  El::LockedView(Z.blocks[block_index], row_offset,
                                 column_offset, Z_block_size, Z_block_size)),
                  Z_times_q(Z_block_size, dx_block_size,
                            Z.blocks[block_index].Grid());
                El::Zero(Z_times_q);
                El::DistMatrix<El::BigFloat> q_Z_q(Z_times_q);
                El::Zero(q_Z_q);

                El::Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL,
                         El::BigFloat(1), Z_sub_block,
                         sdp.bilinear_bases_dist[block_index], El::BigFloat(0),
                         Z_times_q);

                El::Hadamard(Z_times_q, sdp.bilinear_bases_dist[block_index],
                             q_Z_q);

                const size_t dx_row_offset(
                  ((column_block * (column_block + 1)) / 2 + row_block)
                  * dx_block_size);
                El::DistMatrix<El::BigFloat> dx_sub_block(
                  El::View(dx.blocks[jj], dx_row_offset, 0, dx_block_size, 1));

                El::Gemv(El::Orientation::TRANSPOSE, El::BigFloat(-1), q_Z_q,
                         ones, El::BigFloat(1), dx_sub_block);
              }
        }
    }

  // dy = dualObjective - (FreeVarMatrix^T x)
  for(size_t block = 0; block < sdp.free_var_matrix.blocks.size(); ++block)
    {
      El::Gemv(El::OrientationNS::TRANSPOSE, El::BigFloat(-1),
               sdp.free_var_matrix.blocks[block], x.blocks[block],
               El::BigFloat(0), dy.blocks[block]);

      if(block == 0)
        {
          El::Axpy(El::BigFloat(1), sdp.dual_objective_b, dy.blocks[block]);
        }
    }
}
