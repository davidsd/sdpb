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
  //              * Q_Y_Q(column_offset_Y+i,row_offset_Y+j)) / 4

  void multiply_submatrices(const El::DistMatrix<El::BigFloat> &X_inv,
                            const El::DistMatrix<El::BigFloat> &Y,
                            El::DistMatrix<El::BigFloat> &temp,
                            El::DistMatrix<El::BigFloat> &result_submatrix)
  {
    El::Hadamard(X_inv, Y, temp);
    Axpy(El::BigFloat(0.25), temp, result_submatrix);
  }
}

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
                                        schur_complement_block->Grid());

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

                      El::DistMatrix<El::BigFloat> result_sub_matrix(El::View(
                        *schur_complement_block, result_row_offset,
                        result_column_offset, block_size, block_size));
                      Zero(result_sub_matrix);

                      for(size_t parity(0); parity < 2; ++parity)
                        {
                          multiply_submatrices(
                            Q_X_inv_Q[parity][Q_index][column_block_0]
                                     [row_block_1],
                            Q_Y_Q[parity][Q_index][column_block_1][row_block_0],
                            temp, result_sub_matrix);

                          if(block_index == 15)
                            {
                              std::cout.precision(300);
                              std::cout
                                << "result 0: " << block_index << " "
                                << block_size << " " << column_block_0 << " "
                                << row_block_0 << " " << column_block_1 << " "
                                << row_block_1 << " " << parity << " "
                                << El::Norm(result_sub_matrix) << "\n "
                                << El::Norm(
                                     Q_X_inv_Q[parity][Q_index][column_block_0]
                                              [row_block_1])
                                << "\n "
                                << El::Norm(Q_Y_Q[parity][Q_index]
                                                 [column_block_1][row_block_0])
                                << "\n";
                            }

                          multiply_submatrices(
                            Q_X_inv_Q[parity][Q_index][row_block_0]
                                     [row_block_1],
                            Q_Y_Q[parity][Q_index][column_block_1]
                                 [column_block_0],
                            temp, result_sub_matrix);

                          if(block_index == 15)
                            {
                              std::cout.precision(300);
                              std::cout
                                << "result 1: " << block_index << " "
                                << block_size << " " << column_block_0 << " "
                                << row_block_0 << " " << column_block_1 << " "
                                << row_block_1 << " " << parity << " "
                                << El::Norm(result_sub_matrix) << "\n "
                                << El::Norm(
                                     Q_X_inv_Q[parity][Q_index][row_block_0]
                                              [row_block_1])
                                << "\n "
                                << El::Norm(
                                     Q_Y_Q[parity][Q_index][column_block_1]
                                          [column_block_0])
                                << "\n";
                            }

                          if(block_index != 15)
                            {
                              multiply_submatrices(
                                Q_X_inv_Q[parity][Q_index][column_block_0]
                                         [column_block_1],
                                Q_Y_Q[parity][Q_index][row_block_1]
                                     [row_block_0],
                                temp, result_sub_matrix);
                            }
                          if(block_index == 15)
                            {
                              El::Hadamard(
                                Q_X_inv_Q[parity][Q_index][column_block_0]
                                         [column_block_1],
                                Q_Y_Q[parity][Q_index][row_block_1]
                                     [row_block_0],
                                temp);
                              Axpy(El::BigFloat(0.25), temp,
                                   result_sub_matrix);

                              std::cout.precision(300);

                              for(size_t column(0); column < 5; ++column)
                                for(size_t row(0); row < 5; ++row)
                                  {
                                    std::cout
                                      << "temp: " << column << " " << row
                                      << " "
                                      << Q_X_inv_Q[parity][Q_index]
                                                  [column_block_0]
                                                  [column_block_1]
                                                    .Get(column, row)
                                      << "\n "
                                      << Q_Y_Q[parity][Q_index][row_block_1]
                                              [row_block_0]
                                                .Get(column, row)
                                      << "\n " << temp.Get(column, row)
                                      << "\n "
                                      << result_sub_matrix.Get(column, row)
                                      << "\n";
                                  }
                              std::cout
                                << "result 2: " << block_index << " "
                                << block_size << " " << column_block_0 << " "
                                << row_block_0 << " " << column_block_1 << " "
                                << row_block_1 << " " << parity << " "
                                << El::Norm(result_sub_matrix) << "\n "
                                << El::Norm(
                                     Q_X_inv_Q[parity][Q_index][column_block_0]
                                              [column_block_1])
                                << "\n "
                                << El::Norm(Q_Y_Q[parity][Q_index][row_block_1]
                                                 [row_block_0])
                                << "\n " << El::Norm(temp) << "\n";
                            }

                          multiply_submatrices(
                            Q_X_inv_Q[parity][Q_index][row_block_0]
                                     [column_block_1],
                            Q_Y_Q[parity][Q_index][row_block_1][column_block_0],
                            temp, result_sub_matrix);

                          if(block_index == 15)
                            {
                              std::cout.precision(300);
                              std::cout << "result: " << block_index << " "
                                        << block_size << " " << column_block_0
                                        << " " << row_block_0 << " "
                                        << column_block_1 << " " << row_block_1
                                        << " " << parity << " "
                                        << El::Norm(result_sub_matrix) << "\n";
                            }
                        }
                    }
                }
            }
        }

      El::MakeSymmetric(El::UpperOrLower::LOWER, *schur_complement_block);

      if(block_index == 15)
        {
          std::cout.precision(300);
          std::cout << "schur: " << block_index << " " << block_size << " "
                    << dim
                    << " "
                    // // << El::Max(*schur_complement_block) << "\n "
                    // // << El::Min(*schur_complement_block) << "\n "
                    << El::Norm(*schur_complement_block)
                    << "\n "
                    // // << El::SymmetricNorm(El::UpperOrLower::LOWER,
                    // *schur_complement_block) << "\n "
                    // // << El::HermitianNorm(El::UpperOrLower::LOWER,
                    // *schur_complement_block) << "\n "
                    // // << (*schur_complement_block).Get(22,13) << "\n  "
                    // // << (*Q_X_inv_Q_block).at(1).Get(22,13) << "\n  "
                    // // << (*Q_Y_Q_block).at(1).Get(22,13) << " "
                    << El::Norm(Q_X_inv_Q[0][Q_index][0][0]) << "\n "
                    << El::Norm(Q_X_inv_Q[0][Q_index][0][1]) << "\n "
                    << El::Norm(Q_X_inv_Q[0][Q_index][1][1]) << "\n "
                    << El::Norm(Q_Y_Q[0][Q_index][0][0]) << "\n "
                    << El::Norm(Q_Y_Q[0][Q_index][0][1]) << "\n "
                    << El::Norm(Q_Y_Q[0][Q_index][1][1]) << "\n "
                    << "\n";
        }
      ++schur_complement_block;
      ++Q_index;
    }
  schur_complement_timer.stop();
}
