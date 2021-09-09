#include "Zeros.hxx"
#include "../sdp_convert.hxx"

void compute_lambda(const Polynomial_Vector_Matrix &m,
                    const El::Matrix<El::BigFloat> &x,
                    Zeros &zeros)
{
  if(zeros.zeros.empty())
    {
      return;
    }
  El::Matrix<El::BigFloat> interpolation(m.sample_points.size(), zeros.zeros.size());
  for(size_t point_index(0); point_index != m.sample_points.size();
      ++point_index)
    {
      for(size_t zero_index(0); zero_index != zeros.zeros.size(); ++zero_index)
        {
          auto &product(interpolation(point_index, zero_index));
          product = 1;
          for(size_t point_product_index(0);
              point_product_index != m.sample_points.size();
              ++point_product_index)
            {
              if(point_index != point_product_index)
                {
                  product *= (zeros.zeros[zero_index]
                              - m.sample_points[point_product_index])
                             / (m.sample_points[point_index]
                                - m.sample_points[point_product_index]);
                }
            }
        }
    }

  El::Pseudoinverse(interpolation);
  const size_t matrix_block_size(x.Height() / (m.rows * (m.rows + 1) / 2));
  size_t offset(0);
  El::Matrix<El::BigFloat> Lambda(m.rows, m.cols);
  for(int64_t row(0); row != m.rows; ++row)
    for(int64_t column(row); column != m.cols; ++column)
      {
        El::Matrix<El::BigFloat> x_scaled(matrix_block_size, 1);
        for(size_t index(0); index != matrix_block_size; ++index)
          {
            x_scaled(index, 0)
              = x(offset + index, 0) * m.sample_scalings.at(index);
          }
        El::Matrix<El::BigFloat> Lambda_element(
          El::View(Lambda, row, column, 1, 1));
        El::Gemv(El::Orientation::NORMAL,
                 (row == column ? El::BigFloat(1.0) : El::BigFloat(0.5)),
                 interpolation, x_scaled, El::BigFloat(0.0), Lambda_element);
        offset += matrix_block_size;
      }

  El::HermitianEigCtrl<El::BigFloat> hermitian_eig_ctrl;
  /// The default number of iterations is 40.  That is sometimes
  /// not enough, so we bump it up significantly.
  hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.secularCtrl.maxIterations = 400;

  El::Matrix<El::BigFloat> eigenvalues, eigenvectors;
  El::HermitianEig(El::UpperOrLowerNS::UPPER, Lambda, eigenvalues,
                   eigenvectors, hermitian_eig_ctrl);
  // Eigenvalues are sorted, largest at the end
  const size_t num_eigvals(eigenvalues.Height());
  if(eigenvalues(num_eigvals - 1, 0) >= 0)
    {
      zeros.lambda = El::View(eigenvectors, 0, num_eigvals - 1, num_eigvals, 1);
      zeros.lambda *= El::Sqrt(eigenvalues(num_eigvals - 1, 0));
    }
  El::Print(zeros.lambda, "lambda");
  std::cout << "\n";
}
