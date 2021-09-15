#include "Zero.hxx"
#include "../sdp_convert.hxx"

void compute_lambda(const Polynomial_Vector_Matrix &m,
                    const El::Matrix<El::BigFloat> &x,
                    const std::vector<El::BigFloat> &zero_vector,
                    std::vector<Zero> &zeros, El::BigFloat &error)
{
  const size_t matrix_block_size(x.Height() / (m.rows * (m.rows + 1) / 2));

  El::Matrix<El::BigFloat> x_scaled_matrix(matrix_block_size,
                                           (m.rows * (m.rows + 1) / 2));
  size_t row_column(0);
  for(int64_t row(0); row != m.rows; ++row)
    for(int64_t column(row); column != m.cols; ++column)
      {
        for(size_t index(0); index != matrix_block_size; ++index)
          {
            x_scaled_matrix(index, row_column)
              = x(row_column * matrix_block_size + index, 0)
                * m.sample_scalings.at(index);
          }
        ++row_column;
      }
  El::Matrix<El::BigFloat> error_matrix(x_scaled_matrix);

  if(zero_vector.empty())
    {
      error = El::Sqrt(El::Dot(error_matrix, error_matrix));
      return;
    }
  El::Matrix<El::BigFloat> interpolation(m.sample_points.size(),
                                         zero_vector.size());
  for(size_t point_index(0); point_index != m.sample_points.size();
      ++point_index)
    {
      for(size_t zero_index(0); zero_index != zero_vector.size(); ++zero_index)
        {
          auto &product(interpolation(point_index, zero_index));
          product = 1;
          for(size_t point_product_index(0);
              point_product_index != m.sample_points.size();
              ++point_product_index)
            {
              if(point_index != point_product_index)
                {
                  product *= (zero_vector[zero_index]
                              - m.sample_points[point_product_index])
                             / (m.sample_points[point_index]
                                - m.sample_points[point_product_index]);
                }
            }
        }
    }

  El::Matrix<El::BigFloat> roots_fit;
  {
    // This is a copy+paste of El::Pseudoinverse.  We need to set the
    // number of iterations, which is not exposed in the API for
    // El::Pseudoinverse, but is exposed in El::SVD.

    const El::Int m = interpolation.Height();
    const El::Int n = interpolation.Width();
    const El::BigFloat eps = El::limits::Epsilon<El::BigFloat>();

    // Get the SVD
    El::Matrix<El::BigFloat> s;
    El::Matrix<El::BigFloat> U, V;
    El::SVDCtrl<El::BigFloat> ctrl;
    ctrl.overwrite = true;
    ctrl.bidiagSVDCtrl.approach = El::COMPACT_SVD;
    ctrl.bidiagSVDCtrl.tolType = El::RELATIVE_TO_MAX_SING_VAL_TOL;
    ctrl.bidiagSVDCtrl.tol = El::Max(m, n) * eps;
    // The default only does 6 iterations per value.  We need far more
    // because of high precision.
    ctrl.bidiagSVDCtrl.qrCtrl.maxIterPerVal = 100;
    El::SVD(interpolation, U, s, V, ctrl);

    // Scale U with the inverted (nonzero) singular values, U := U / Sigma
    El::DiagonalSolve(El::RIGHT, El::NORMAL, s, U);

    // Form pinvA = (U Sigma V^H)^H = V (U Sigma)^H
    El::Gemm(El::NORMAL, El::ADJOINT, El::BigFloat(1), V, U, roots_fit);
  }

  for(size_t zero_index(0); zero_index != zero_vector.size(); ++zero_index)
    {
      El::Matrix<El::BigFloat> Lambda(m.rows, m.cols);
      El::Matrix<El::BigFloat> fit(roots_fit.Height(), 1);
      size_t row_column(0);
      for(int64_t row(0); row != m.rows; ++row)
        for(int64_t column(row); column != m.cols; ++column)
          {
            El::Matrix<El::BigFloat> roots_view(
              El::View(roots_fit, zero_index, 0, 1, roots_fit.Width()));
            El::Matrix<El::BigFloat> Lambda_view(
              El::View(Lambda, row, column, 1, 1));

            El::Gemv(
              El::Orientation::NORMAL,
              (row == column ? El::BigFloat(1.0) : El::BigFloat(0.5)),
              roots_view,
              El::View(x_scaled_matrix, 0, row_column, matrix_block_size, 1),
              El::BigFloat(0.0), Lambda_view);
            Lambda(column, row) = Lambda(row, column);
            ++row_column;
          }

      El::HermitianEigCtrl<El::BigFloat> hermitian_eig_ctrl;
      /// The default number of iterations is 40.  That is sometimes
      /// not enough, so we bump it up significantly.
      hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.secularCtrl.maxIterations = 400;

      El::Matrix<El::BigFloat> eigenvalues, eigenvectors;
      El::HermitianEig(El::UpperOrLowerNS::UPPER, Lambda, eigenvalues,
                       eigenvectors, hermitian_eig_ctrl);
      // Eigenvalues are sorted, largest at the end.  Only add a zero
      // if max_eigenvalue > 0.
      const size_t num_eigvals(eigenvalues.Height());
      if(eigenvalues(num_eigvals - 1, 0) >= 0)
        {
          zeros.emplace_back(zero_vector[zero_index]);
          zeros.back().lambda
            = El::View(eigenvectors, 0, num_eigvals - 1, num_eigvals, 1);
          zeros.back().lambda *= El::Sqrt(eigenvalues(num_eigvals - 1, 0));
        }
    }
  error = El::Sqrt(El::Dot(error_matrix, error_matrix));
}
