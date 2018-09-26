#include "../../../Block_Diagonal_Matrix.hxx"
#include "../../../../Block_Info.hxx"
#include "../../../../Timers.hxx"

// Compute the SchurComplement matrix using BilinearPairingsXInv and
// BilinearPairingsY and the formula
//
//   S_{(j,r1,s1,k1), (j,r2,s2,k2)} = \sum_{b \in blocks[j]}
//          (1/4) (BilinearPairingsXInv_{ej s1 + k1, ej r2 + k2}*
//                 BilinearPairingsY_{ej s2 + k2, ej r1 + k1} +
//                 swaps (r1 <-> s1) and (r2 <-> s2))
//
// where ej = d_j + 1.

namespace
{
  // Elementwise multiplication of two submatrices.
  //
  // result(i,j)=bilinear_X_inv(column_offset_X+i,row_offset_X+j)
  //             + bilinear_Y_inv(column_offset_Y+i,row_offset_Y+j)

  void
  multiply_submatrices(const El::DistMatrix<El::BigFloat> &bilinear_X_inv,
                       const El::DistMatrix<El::BigFloat> &bilinear_Y,
                       const size_t &block_size, const size_t &column_offset_X,
                       const size_t &row_offset_X,
                       const size_t &column_offset_Y,
                       const size_t &row_offset_Y,
                       El::DistMatrix<El::BigFloat> &temp,
                       El::DistMatrix<El::BigFloat> &result_submatrix)
  {
    El::DistMatrix<El::BigFloat> X_submatrix(El::LockedView(
      bilinear_X_inv, column_offset_X, row_offset_X, block_size, block_size)),
      Y_submatrix(El::LockedView(bilinear_Y, column_offset_Y, row_offset_Y,
                                 block_size, block_size)),
      Y_transpose(block_size, block_size, bilinear_Y.Grid());

    // The overall result is symmetric, but this particular block may
    // not be because it is offset from the diagonal.

    // It does not help to use TransposeAxpy().  Hadamard() sets the
    // result, it does not add it in.
    El::Transpose(Y_submatrix, Y_transpose);
    if(X_submatrix.ColAlign() == Y_transpose.ColAlign()
       && X_submatrix.RowAlign() == Y_transpose.RowAlign())
      {
        El::Hadamard(X_submatrix, Y_transpose, temp);
      }
    else
      {
        El::DistMatrix<El::BigFloat> Y_transpose_aligned(
          block_size, block_size, Y_transpose.Grid());
        Y_transpose_aligned.AlignWith(X_submatrix);
        El::Copy(Y_transpose, Y_transpose_aligned);
        El::Hadamard(X_submatrix, Y_transpose_aligned, temp);
      }
    Axpy(El::BigFloat(0.25), temp, result_submatrix);
  }
}

void compute_schur_complement(
  const Block_Info &block_info,
  const Block_Diagonal_Matrix &bilinear_pairings_X_inv,
  const Block_Diagonal_Matrix &bilinear_pairings_Y,
  Block_Diagonal_Matrix &schur_complement, Timers &timers)
{
  auto &schur_complement_timer(timers.add_and_start(
    "run.step.initializeSchurComplementSolver.schur_complement"));

  auto schur_complement_block(schur_complement.blocks.begin());
  auto bilinear_pairings_X_inv_block(bilinear_pairings_X_inv.blocks.begin());
  auto bilinear_pairings_Y_block(bilinear_pairings_Y.blocks.begin());
  for(auto &block_index : block_info.block_indices)
    {
      const size_t block_size(block_info.degrees[block_index] + 1);

      for(size_t column_block_0 = 0;
          column_block_0 < block_info.dimensions[block_index];
          ++column_block_0)
        {
          const size_t column_offset_0(column_block_0 * block_size);
          for(size_t row_block_0 = 0; row_block_0 <= column_block_0;
              ++row_block_0)
            {
              const size_t row_offset_0(row_block_0 * block_size);
              size_t result_row_offset(
                ((column_block_0 * (column_block_0 + 1)) / 2 + row_block_0)
                * block_size);

              for(size_t column_block_1 = 0;
                  column_block_1 < block_info.dimensions[block_index];
                  ++column_block_1)
                {
                  const size_t column_offset_1(column_block_1 * block_size);
                  for(size_t row_block_1 = 0; row_block_1 <= column_block_1;
                      ++row_block_1)
                    {
                      const size_t row_offset_1(row_block_1 * block_size);
                      size_t result_column_offset(
                        ((column_block_1 * (column_block_1 + 1)) / 2
                         + row_block_1)
                        * block_size);

                      El::DistMatrix<El::BigFloat> result_sub_matrix(El::View(
                        *schur_complement_block, result_row_offset,
                        result_column_offset, block_size, block_size));
                      Zero(result_sub_matrix);

                      El::DistMatrix<El::BigFloat> temp(
                        block_size, block_size, result_sub_matrix.Grid());

                      auto bilinear_X_inv(bilinear_pairings_X_inv_block);
                      auto bilinear_Y(bilinear_pairings_Y_block);
                      for(size_t bb(2 * block_index); bb < 2 * block_index + 2;
                          ++bb)
                        {
                          multiply_submatrices(
                            *bilinear_X_inv, *bilinear_Y, block_size,
                            column_offset_0, row_offset_1, column_offset_1,
                            row_offset_0, temp, result_sub_matrix);
                          multiply_submatrices(
                            *bilinear_X_inv, *bilinear_Y, block_size,
                            row_offset_0, row_offset_1, column_offset_1,
                            column_offset_0, temp, result_sub_matrix);
                          multiply_submatrices(
                            *bilinear_X_inv, *bilinear_Y, block_size,
                            column_offset_0, column_offset_1, row_offset_1,
                            row_offset_0, temp, result_sub_matrix);
                          multiply_submatrices(
                            *bilinear_X_inv, *bilinear_Y, block_size,
                            row_offset_0, column_offset_1, row_offset_1,
                            column_offset_0, temp, result_sub_matrix);
                          ++bilinear_X_inv;
                          ++bilinear_Y;
                        }
                    }
                }
            }
        }

      El::MakeSymmetric(El::UpperOrLower::LOWER, *schur_complement_block);
      ++schur_complement_block;
      ++bilinear_pairings_X_inv_block;
      ++bilinear_pairings_X_inv_block;
      ++bilinear_pairings_Y_block;
      ++bilinear_pairings_Y_block;
    }
  schur_complement_timer.stop();
}
