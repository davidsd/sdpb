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
    std::vector<std::vector<std::vector<El::Matrix<El::BigFloat>>>>, 2>
    &Q_X_inv_Q,
  const std::array<
    std::vector<std::vector<std::vector<El::Matrix<El::BigFloat>>>>, 2> &Q_Y_Q,
  Block_Diagonal_Matrix &schur_complement, Timers &timers)
{
  auto &schur_complement_timer(timers.add_and_start(
    "run.step.initializeSchurComplementSolver.schur_complement"));

  auto schur_complement_block(schur_complement.blocks.begin());
  size_t Q_index(0);
  // Put these BigFloats at the beginning to avoid memory churn
  El::BigFloat element, product;
  for(auto &block_index : block_info.block_indices)
    {
      const size_t block_size(block_info.num_points[block_index]),
        dim(block_info.dimensions[block_index]);

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
                      size_t result_column_offset(
                        ((column_block_1 * (column_block_1 + 1)) / 2
                         + row_block_1)
                        * block_size);

                      for(size_t row(0); row < block_size; ++row)
                        {
                          for(size_t column(0); column < block_size; ++column)
                            {
                              size_t global_row(row + result_row_offset),
                                global_column(column + result_column_offset);
                              if(schur_complement_block->IsLocal(global_row, global_column))
                                {
                                  element.Zero();
                                  for(size_t parity(0); parity < 2; ++parity)
                                    {
                                      // Do this the hard way to avoid memory
                                      // allocations
                                      product
                                        = Q_X_inv_Q[parity][Q_index]
                                                   [column_block_0][row_block_1]
                                                     .CRef(row, column);
                                      product
                                        *= Q_Y_Q[parity][Q_index]
                                                [column_block_1][row_block_0]
                                                  .CRef(row, column);
                                      element += product;

                                      product
                                        = Q_X_inv_Q[parity][Q_index]
                                                   [row_block_0][row_block_1]
                                                     .CRef(row, column);
                                      product
                                        *= Q_Y_Q[parity][Q_index]
                                                [column_block_1][column_block_0]
                                                  .CRef(row, column);
                                      element += product;

                                      product = Q_X_inv_Q[parity][Q_index]
                                                         [column_block_0]
                                                         [column_block_1]
                                                           .CRef(row, column);
                                      product
                                        *= Q_Y_Q[parity][Q_index][row_block_1]
                                                [row_block_0]
                                                  .CRef(row, column);
                                      element += product;

                                      product
                                        = Q_X_inv_Q[parity][Q_index]
                                                   [row_block_0][column_block_1]
                                                     .CRef(row, column);
                                      product
                                        *= Q_Y_Q[parity][Q_index][row_block_1]
                                                [column_block_0]
                                                  .CRef(row, column);
                                      element += product;
                                    }
                                  element /= 4;

                                  int64_t local_row(
                                    schur_complement_block->LocalRow(
                                      row + result_row_offset)),
                                    local_column(
                                      schur_complement_block->LocalCol(
                                        column + result_column_offset));
                                  schur_complement_block->SetLocal(
                                    local_row, local_column, element);
                                }
                            }
                        }
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
