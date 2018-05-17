#include "../../../../SDP_Solver.hxx"

// Compute the vectors r_x and r_y on the right-hand side of the Schur
// complement equation:
//
// {{S, -B}, {B^T, 0}} . {dx, dy} = {r_x, r_y}
//
// where S = SchurComplement and B = FreeVarMatrix.  Specifically,
//
// r_x[p] = -dual_residues[p] - Tr(A_p Z)              for 0 <= p < P
// r_y[n] = dualObjective[n] - (FreeVarMatrix^T x)_n  for 0 <= n < N
//
// where P = primalObjective.size(), N = dualObjective.size()
//
// Inputs:
// - sdp, an SDP
// - dual_residues, a Vector of length P
// - Z = X^{-1} (PrimalResidues Y - R)
// - x, a vector of length P
// Outputs:
// - r_x, a Vector of length P
// - r_y, a Vector of length N
//

Real bilinear_block_pairing(const Real *v, const int dim, const Matrix &A,
                            const int blockRow, const int blockCol);

void compute_schur_RHS(const SDP &sdp,
                       const Block_Vector &dual_residues_elemental,
                       const Block_Diagonal_Matrix &Z,
                       const Block_Vector &x_elemental,
                       Block_Vector &r_x_elemental,
                       El::DistMatrix<El::BigFloat> &r_y_elemental)
{
  for(size_t jj = 0; jj < sdp.dimensions.size(); ++jj)
    {
      // r_x = -dual_residues
      r_x_elemental.blocks[jj] = dual_residues_elemental.blocks[jj];
      r_x_elemental.blocks[jj] *= -1;
      const size_t r_x_block_size(sdp.degrees[jj] + 1);

      // r_x[p] -= Tr(A_p Z)
      // Not sure whether it is better to first loop over blocks in
      // the result or over sub-blocks in Z
      for(size_t block_index = 2 * jj; block_index < 2 * jj + 2; ++block_index)
        {
          const size_t Z_block_size(
            sdp.bilinear_bases_elemental_dist[block_index].Height());
          El::DistMatrix<El::BigFloat> ones;
          El::Ones(ones, Z_block_size, 1);

          for(size_t column_block = 0; column_block < sdp.dimensions[jj];
              ++column_block)
            for(size_t row_block = 0; row_block <= column_block; ++row_block)
              {
                size_t column_offset(column_block * Z_block_size),
                  row_offset(row_block * Z_block_size);

                El::DistMatrix<El::BigFloat> Z_sub_block(
                  El::LockedView(Z.blocks_elemental[block_index], row_offset,
                                 column_offset, Z_block_size, Z_block_size)),
                  Z_times_q(Z_block_size, r_x_block_size),
                  q_Z_q(Z_block_size, r_x_block_size);
                El::Zero(Z_times_q);
                El::Zero(q_Z_q);

                El::Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL,
                         El::BigFloat(1), Z_sub_block,
                         sdp.bilinear_bases_elemental_dist[block_index],
                         El::BigFloat(0), Z_times_q);

                El::Hadamard(Z_times_q,
                             sdp.bilinear_bases_elemental_dist[block_index],
                             q_Z_q);

                const size_t r_x_row_offset(
                  ((column_block * (column_block + 1)) / 2 + row_block)
                  * r_x_block_size);
                El::DistMatrix<El::BigFloat> r_x_sub_block(
                  El::View(r_x_elemental.blocks[jj], r_x_row_offset, 0,
                           r_x_block_size, 1));

                El::Gemv(El::Orientation::TRANSPOSE, El::BigFloat(-1), q_Z_q,
                         ones, El::BigFloat(1), r_x_sub_block);
              }
        }
    }

  // r_y = dualObjective - (FreeVarMatrix^T x)
  r_y_elemental = sdp.dual_objective_b_elemental;
  for(size_t b = 0; b < sdp.free_var_matrix_elemental.blocks.size(); ++b)
    {
      El::Gemv(El::OrientationNS::TRANSPOSE, El::BigFloat(-1),
               sdp.free_var_matrix_elemental.blocks[b], x_elemental.blocks[b],
               El::BigFloat(1), r_y_elemental);
    }
}
