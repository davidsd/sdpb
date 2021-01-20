#include "../../../../Block_Diagonal_Matrix.hxx"
#include "../../../../Block_Info.hxx"
#include "../../../../../Timers.hxx"

// Compute the SchurComplement matrix using Q_X_inv_Q and
// Q_Y_Q and the formula
//
//   S_{(j,r1,s1,k1), (j,r2,s2,k2)} = \sum_{b \in blocks[j]}
//          (1/4) (Q_X_inv_Q_{ej s1 + k1, ej r2 + k2}*
//                 Q_Y_Q_{ej s2 + k2, ej r1 + k1} +
//                 swaps (r1 <-> s1) and (r2 <-> s2))
//
// where ej = d_j + 1.

namespace
{
  // Elementwise multiplication of two submatrices.
  //
  // result(i,j)=(Q_X_inv_Q(column_offset_X+i,row_offset_X+j)
  //              * Q_Y_Q(column_offset_Y+i,row_offset_Y+j))/4

  void
  multiply_submatrices(const El::DistMatrix<El::BigFloat> &X_inv,
                       const El::DistMatrix<El::BigFloat> &Y,
                       const size_t &block_size, const size_t &column_offset_X,
                       const size_t &row_offset_X,
                       const size_t &column_offset_Y,
                       const size_t &row_offset_Y,
                       El::DistMatrix<El::BigFloat> &temp,
                       El::DistMatrix<El::BigFloat> &result_submatrix)
  {
    const El::DistMatrix<El::BigFloat> X_submatrix(El::LockedView(
      X_inv, column_offset_X, row_offset_X, block_size, block_size)),
      Y_submatrix(El::LockedView(Y, column_offset_Y, row_offset_Y, block_size,
                                 block_size));
    if(X_submatrix.ColAlign() == Y_submatrix.ColAlign()
       && X_submatrix.RowAlign() == Y_submatrix.RowAlign())
      {
        El::Hadamard(X_submatrix, Y_submatrix, temp);
      }
    else
      {
        El::DistMatrix<El::BigFloat> aligned(
          block_size, block_size, Y_submatrix.Grid());
        aligned.AlignWith(X_submatrix);
        El::Copy(Y_submatrix, aligned);
        El::Hadamard(X_submatrix, aligned, temp);
      }
    Axpy(El::BigFloat(0.25), temp, result_submatrix);
  }
}

void compute_schur_complement(const Block_Info &block_info,
                              const Block_Diagonal_Matrix &Q_X_inv_Q,
                              const Block_Diagonal_Matrix &Q_Y_Q,
                              Block_Diagonal_Matrix &schur_complement,
                              Timers &timers)
{
  auto &schur_complement_timer(timers.add_and_start(
    "run.step.initializeSchurComplementSolver.schur_complement"));

  auto schur_complement_block(schur_complement.blocks.begin());
  auto Q_X_inv_Q_block(Q_X_inv_Q.blocks.begin());
  auto Q_Y_Q_block(Q_Y_Q.blocks.begin());
  for(auto &block_index : block_info.block_indices)
    {
      const size_t block_size(block_info.num_points[block_index]);

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

                      auto X_inv(Q_X_inv_Q_block);
                      auto Y(Q_Y_Q_block);
                      for(size_t parity(0); parity < 2; ++parity)
                        {
                          multiply_submatrices(*X_inv, *Y, block_size,
                                               column_offset_0, row_offset_1,
                                               column_offset_1, row_offset_0,
                                               temp, result_sub_matrix);
                          multiply_submatrices(
                            *X_inv, *Y, block_size, row_offset_0, row_offset_1,
                            column_offset_1, column_offset_0, temp,
                            result_sub_matrix);
                          multiply_submatrices(
                            *X_inv, *Y, block_size, column_offset_0,
                            column_offset_1, row_offset_1, row_offset_0, temp,
                            result_sub_matrix);
                          multiply_submatrices(*X_inv, *Y, block_size,
                                               row_offset_0, column_offset_1,
                                               row_offset_1, column_offset_0,
                                               temp, result_sub_matrix);
                          ++X_inv;
                          ++Y;
                        }
                    }
                }
            }
        }

      El::MakeSymmetric(El::UpperOrLower::LOWER, *schur_complement_block);
      ++schur_complement_block;
      ++Q_X_inv_Q_block;
      ++Q_X_inv_Q_block;
      ++Q_Y_Q_block;
      ++Q_Y_Q_block;
    }
  schur_complement_timer.stop();
}
