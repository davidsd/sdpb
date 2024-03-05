#include "sdp_solve/Block_Diagonal_Matrix.hxx"
#include "sdp_solve/Block_Info.hxx"
#include "sdpb_util/Timers/Timers.hxx"

// Compute the SchurComplement matrix using A_X_inv and
// A_Y and the formula
//
//   S_{(j,r1,s1,k1), (j,r2,s2,k2)} = \sum_{b \in blocks[j]}
//          (1/4) (A_X_inv_{ej s1 + k1, ej r2 + k2}*
//                 A_Y_{ej s2 + k2, ej r1 + k1} +
//                 swaps (r1 <-> s1) and (r2 <-> s2))
//
// where ej = d_j + 1.

void compute_schur_complement(
  const Block_Info &block_info,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_X_inv,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
  Block_Diagonal_Matrix &schur_complement, Timers &timers)
{
  Scoped_Timer schur_complement_timer(timers, "schur_complement");

  auto schur_complement_block(schur_complement.blocks.begin());
  size_t Q_index(0);
  // Put these BigFloats at the beginning to avoid memory churn
  El::BigFloat element, product;
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
                                    = A_X_inv[parity][Q_index][column_block_0]
                                             [row_block_1]
                                               .GetLocalCRef(row, column);
                                  product *= A_Y[parity][Q_index]
                                                [column_block_1][row_block_0]
                                                  .GetLocalCRef(row, column);
                                  element += product;

                                  product
                                    = A_X_inv[parity][Q_index][row_block_0]
                                             [row_block_1]
                                               .GetLocalCRef(row, column);
                                  product
                                    *= A_Y[parity][Q_index][column_block_1]
                                          [column_block_0]
                                            .GetLocalCRef(row, column);
                                  element += product;

                                  product
                                    = A_X_inv[parity][Q_index][column_block_0]
                                             [column_block_1]
                                               .GetLocalCRef(row, column);
                                  product *= A_Y[parity][Q_index][row_block_1]
                                                [row_block_0]
                                                  .GetLocalCRef(row, column);
                                  element += product;

                                  product
                                    = A_X_inv[parity][Q_index][row_block_0]
                                             [column_block_1]
                                               .GetLocalCRef(row, column);
                                  product *= A_Y[parity][Q_index][row_block_1]
                                                [column_block_0]
                                                  .GetLocalCRef(row, column);
                                  element += product;
                                }
                              element /= 4;
                              temp_result.SetLocal(row, column, element);
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
}
