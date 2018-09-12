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
void compute_dual_residues(const Block_Info &block_info, const SDP &sdp,
                           const Block_Vector &y,
                           const Block_Diagonal_Matrix &bilinear_pairings_Y,
                           Block_Vector &dual_residues)
{
  auto dual_residues_block(dual_residues.blocks.begin());
  auto primal_objective_c_block(sdp.primal_objective_c.blocks.begin());
  auto y_block(y.blocks.begin());
  auto free_var_matrix_block(sdp.free_var_matrix.blocks.begin());
  auto bilinear_pairings_Y_block(bilinear_pairings_Y.blocks.begin());
  // dualResidues[p] = -Tr(A_p Y)
  for(auto &block_index : block_info.block_indices)
    {
      Zero(*dual_residues_block);
      const size_t block_size(block_info.degrees[block_index] + 1);

      // Not sure whether it is better to first loop over blocks in
      // the result or over sub-blocks in bilinear_pairings_Y
      for(size_t bb = 2 * block_index; bb < 2 * block_index + 2; ++bb)
        {
          for(size_t column_block = 0;
              column_block < block_info.dimensions[block_index];
              ++column_block)
            for(size_t row_block = 0; row_block <= column_block; ++row_block)
              {
                size_t column_offset(column_block * block_size),
                  row_offset(row_block * block_size);

                El::DistMatrix<El::BigFloat> lower_diagonal(El::GetDiagonal(
                  El::LockedView(*bilinear_pairings_Y_block, row_offset,
                                 column_offset, block_size, block_size)));

                size_t residue_row_offset(
                  ((column_block * (column_block + 1)) / 2 + row_block)
                  * block_size);

                El::DistMatrix<El::BigFloat> residue_sub_block(El::View(
                  *dual_residues_block, residue_row_offset, 0, block_size, 1));

                // Do a little less work when on a diagonal block.
                if(column_offset == row_offset)
                  {
                    El::Axpy(El::BigFloat(-1.0), lower_diagonal,
                             residue_sub_block);
                  }
                else
                  {
                    El::Axpy(El::BigFloat(-0.5), lower_diagonal,
                             residue_sub_block);

                    El::DistMatrix<El::BigFloat> upper_diagonal(
                      El::GetDiagonal(El::LockedView(
                        *bilinear_pairings_Y_block, column_offset, row_offset,
                        block_size, block_size)));

                    El::Axpy(El::BigFloat(-0.5), upper_diagonal,
                             residue_sub_block);
                  }
              }
          ++bilinear_pairings_Y_block;
        }
      // dualResidues -= FreeVarMatrix * y
      Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL, El::BigFloat(-1),
           *free_var_matrix_block, *y_block, El::BigFloat(1),
           *dual_residues_block);
      // dualResidues += primalObjective
      Axpy(El::BigFloat(1), *primal_objective_c_block, *dual_residues_block);

      ++primal_objective_c_block;
      ++y_block;
      ++free_var_matrix_block;
      ++dual_residues_block;
    }
}
