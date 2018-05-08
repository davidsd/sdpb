#include "../../SDP_Solver.hxx"

// dualResidues[p] = primalObjective[p] - Tr(A_p Y) - (FreeVarMatrix y)_p,
// for 0 <= p < primalObjective.size()
//
// The pairings Tr(A_p Y) can be written in terms of BilinearPairingsY:
//
//   Tr(A_(j,r,s,k) Y) = \sum_{b \in blocks[j]}
//                       (1/2) (BilinearPairingsY_{ej r + k, ej s + k} +
//                              swap (r <-> s))
// where ej = d_j + 1.
//
// Inputs: sdp, y, BilinearPairingsY
// Output: dualResidues (overwriten)
//
void compute_dual_residues(const SDP &sdp, const Vector &y,
                           const Block_Diagonal_Matrix &bilinear_pairings_Y,
                           Vector &dual_residues)
{
  for(size_t j = 0; j < sdp.dimensions.size(); j++)
    {
      const int ej = sdp.degrees[j] + 1;

      for(auto &t : sdp.constraint_indices[j])
        {
          const int p = t.p;
          const int ej_r = t.r * ej;
          const int ej_s = t.s * ej;
          const int k = t.k;

          // dualResidues[p] = -Tr(A_p Y)
          dual_residues[p] = 0;
          for(auto &b : sdp.blocks[j])
            {
              dual_residues[p]
                -= bilinear_pairings_Y.blocks[b].elt(ej_r + k, ej_s + k);
              dual_residues[p]
                -= bilinear_pairings_Y.blocks[b].elt(ej_s + k, ej_r + k);
            }
          dual_residues[p] /= 2;

          // dualResidues[p] = -Tr(A_p Y) - (FreeVarMatrix y)_p
          for(size_t n = 0; n < sdp.free_var_matrix.cols; n++)
            {
              dual_residues[p] -= sdp.free_var_matrix.elt(p, n) * y[n];
            }

          // dualResidues[p] = primalObjective[p] - Tr(A_p Y) - (FreeVarMatrix
          // y)_p
          dual_residues[p] += sdp.primal_objective_c[p];
        }
    }
}

void compute_dual_residues(const SDP &sdp,
                           const El::DistMatrix<El::BigFloat> &y,
                           const Block_Diagonal_Matrix &bilinear_pairings_Y,
                           Block_Matrix &dual_residues)
{
  // dualResidues[p] = -Tr(A_p Y)
  for(size_t jj = 0; jj < sdp.dimensions.size(); ++jj)
    {
      Zero(dual_residues.blocks[jj]);
      const size_t block_size(sdp.degrees[jj] + 1);

      // Not sure whether it is better to first loop over blocks in
      // the result or over sub-blocks in bilinear_pairings_Y
      for(size_t block_index = 2 * jj; block_index < 2 * jj + 2; ++block_index)
        {
          for(size_t column_block = 0; column_block < sdp.dimensions[jj];
              ++column_block)
            for(size_t row_block = 0; row_block <= column_block; ++row_block)
              {
                size_t column_offset(column_block * block_size),
                  row_offset(row_block * block_size);

                El::DistMatrix<El::BigFloat> lower_diagonal(
                  El::GetDiagonal(El::LockedView(
                    bilinear_pairings_Y.blocks_elemental[block_index],
                    row_offset, column_offset, block_size, block_size))),
                  upper_diagonal(El::GetDiagonal(El::LockedView(
                    bilinear_pairings_Y.blocks_elemental[block_index],
                    column_offset, row_offset, block_size, block_size)));

                size_t residue_row_offset(
                  ((column_block * (column_block + 1)) / 2 + row_block)
                  * block_size);

                El::DistMatrix<El::BigFloat> residue_sub_block(
                  El::View(dual_residues.blocks[jj], residue_row_offset, 0,
                           block_size, 1));

                // FIXME: This does twice as much work as needed when
                // column_block==row_block.
                El::Axpy(El::BigFloat(-0.5), lower_diagonal,
                         residue_sub_block);
                El::Axpy(El::BigFloat(-0.5), upper_diagonal,
                         residue_sub_block);
              }
        }
    }

  for(size_t b = 0; b < dual_residues.blocks.size(); ++b)
    {
      // dualResidues -= FreeVarMatrix * y
      Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL, El::BigFloat(-1),
           sdp.free_var_matrix_elemental.blocks[b], y, El::BigFloat(1),
           dual_residues.blocks[b]);
      // dualResidues += primalObjective
      Axpy(El::BigFloat(1), sdp.primal_objective_c_elemental.blocks[b],
           dual_residues.blocks[b]);
    }
}
