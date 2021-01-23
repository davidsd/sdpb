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

void compute_schur_complement(
  const Block_Info &block_info,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &Q_X_inv_Q,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &Q_Y_Q,
  Block_Diagonal_Matrix &schur_complement, Timers &timers)
{
  auto &schur_complement_timer(timers.add_and_start(
    "run.step.initializeSchurComplementSolver.schur_complement"));

  auto schur_complement_block(schur_complement.blocks.begin());
  size_t Q_index(0);
  for(auto &block_index : block_info.block_indices)
    {
      const size_t block_size(block_info.num_points[block_index]),
        dim(block_info.dimensions[block_index]);

      El::DistMatrix<El::BigFloat> temp(block_size, block_size,
                                        schur_complement_block->Grid()),
        temp_result(temp);

      for(size_t column_block_0 = 0; column_block_0 < dim; ++column_block_0)
        {
          for(size_t row_block_0 = 0; row_block_0 <= column_block_0;
              ++row_block_0)
            {
              size_t result_row_offset(
                ((column_block_0 * (column_block_0 + 1)) / 2 + row_block_0)
                * block_size);

              for(size_t column_block_1 = 0; column_block_1 < dim;
                  ++column_block_1)
                {
                  for(size_t row_block_1 = 0; row_block_1 <= column_block_1;
                      ++row_block_1)
                    {
                      El::BigFloat element, product;
                      for(int64_t row(0); row < temp_result.LocalHeight();
                          ++row)
                        {
                          for(int64_t column(0);
                              column < temp_result.LocalWidth(); ++column)
                            {
                              element.Zero();
                              for(size_t parity(0); parity < 2; ++parity)
                                {
                                  // Do this the hard way to avoid memory
                                  // allocations
                                  product
                                    = Q_X_inv_Q[parity][Q_index]
                                               [column_block_0][row_block_1]
                                                 .GetLocal(row, column);
                                  product *= Q_Y_Q[parity][Q_index]
                                                  [column_block_1][row_block_0]
                                                    .GetLocal(row, column);
                                  element += product;

                                  product = Q_X_inv_Q[parity][Q_index]
                                                     [row_block_0][row_block_1]
                                                       .GetLocal(row, column);
                                  product
                                    *= Q_Y_Q[parity][Q_index][column_block_1]
                                            [column_block_0]
                                              .GetLocal(row, column);
                                  element += product;

                                  product
                                    = Q_X_inv_Q[parity][Q_index]
                                               [column_block_0][column_block_1]
                                                 .GetLocal(row, column);
                                  product *= Q_Y_Q[parity][Q_index]
                                                  [row_block_1][row_block_0]
                                                    .GetLocal(row, column);
                                  element += product;

                                  product
                                    = Q_X_inv_Q[parity][Q_index][row_block_0]
                                               [column_block_1]
                                                 .GetLocal(row, column);
                                  product *= Q_Y_Q[parity][Q_index]
                                                  [row_block_1][column_block_0]
                                                    .GetLocal(row, column);
                                  element += product;
                                }
                              element /= 4;
                              temp_result.SetLocal(column, row, element);
                            }
                        }
                      size_t result_column_offset(
                        ((column_block_1 * (column_block_1 + 1)) / 2
                         + row_block_1)
                        * block_size);

                      El::DistMatrix<El::BigFloat> result_submatrix(El::View(
                        *schur_complement_block, result_row_offset,
                        result_column_offset, block_size, block_size));

                      El::Copy(temp_result, result_submatrix);
                    }
                }
            }
        }

      El::MakeSymmetric(El::UpperOrLower::LOWER, *schur_complement_block);
      ++schur_complement_block;
      ++Q_index;
    }
  schur_complement_timer.stop();
}
