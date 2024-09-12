#pragma once

#include "pmp/PMP_Info.hxx"
#include "spectrum/Zero.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/Damped_Rational.hxx"
#include "sdpb_util/Timers/Timers.hxx"

inline void
compute_lambda(const PVM_Info &pvm_info, const El::Matrix<El::BigFloat> &x,
               const std::vector<El::BigFloat> &zero_values,
               Zeros &spectrum_block, Timers &timers)
{
  Scoped_Timer timer(timers, "compute_lambda");

  const auto num_rows = pvm_info.dim;
  const auto &sample_points = pvm_info.sample_points;

  auto &zeros = spectrum_block.zeros;
  auto &error = spectrum_block.error;

  const size_t matrix_block_size(x.Height() / (num_rows * (num_rows + 1) / 2));

  // U_{j,k} in Eq. (A.11)
  // hereafter equation numbers are from https://arxiv.org/pdf/1612.08471
  El::Matrix<El::BigFloat> x_scaled_matrix(matrix_block_size,
                                           (num_rows * (num_rows + 1) / 2));
  {
    size_t row_column(0);
    for(size_t column(0); column != num_rows; ++column)
      for(size_t row(0); row <= column; ++row)
        {
          for(size_t index(0); index != matrix_block_size; ++index)
            {
              x_scaled_matrix(index, row_column)
                = x(row_column * matrix_block_size + index, 0)
                  * pvm_info.reduced_sample_scalings.at(index);
            }
          ++row_column;
        }
  }
  El::Matrix<El::BigFloat> error_matrix(x_scaled_matrix);

  if(zero_values.empty())
    {
      error = El::Sqrt(El::Dot(error_matrix, error_matrix));
      return;
    }

  // Lagrange interpolation coefficients L(\tau, x_k^{(j)}, Eq. (A.15)
  El::Matrix<El::BigFloat> interpolation(sample_points.size(),
                                         zero_values.size());
  for(size_t point_index(0); point_index != sample_points.size();
      ++point_index)
    {
      for(size_t zero_index(0); zero_index != zero_values.size(); ++zero_index)
        {
          auto &product(interpolation(point_index, zero_index));
          product = 1;
          for(size_t point_product_index(0);
              point_product_index != sample_points.size();
              ++point_product_index)
            {
              if(point_index != point_product_index)
                {
                  product *= (zero_values[zero_index]
                              - sample_points[point_product_index])
                             / (sample_points[point_index]
                                - sample_points[point_product_index]);
                }
            }
        }
    }

  // Solve Eq. (A.15) for V using least-squares fit.
  // roots_fit is L^{-1}
  El::Matrix<El::BigFloat> roots_fit;
  {
    // This is mostly a copy+paste of El::Pseudoinverse.  We need to
    // set the number of iterations, which is not exposed in the API
    // for El::Pseudoinverse, but is exposed in El::SVD.

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
    El::Matrix<El::BigFloat> interpolation_copy(interpolation);
    El::SVD(interpolation_copy, U, s, V, ctrl);

    // Scale U with the inverted (nonzero) singular values, U := U / Sigma
    El::DiagonalSolve(El::RIGHT, El::NORMAL, s, U);

    // Form pinvA = (U Sigma V^H)^H = V (U Sigma)^H
    El::Gemm(El::NORMAL, El::ADJOINT, El::BigFloat(1), V, U, roots_fit);
  }

  for(size_t zero_index(0); zero_index != zero_values.size(); ++zero_index)
    {
      //  V_{j,\tau} from Eq. (A.15)
      El::Matrix<El::BigFloat> Lambda(num_rows, num_rows);
      El::Matrix<El::BigFloat> fit(roots_fit.Height(), 1);
      size_t row_column(0);
      for(size_t column(0); column != num_rows; ++column)
        for(size_t row(0); row <= column; ++row)
          {
            El::Matrix<El::BigFloat> roots_view(
              El::View(roots_fit, zero_index, 0, 1, roots_fit.Width()));
            El::Matrix<El::BigFloat> Lambda_view(
              El::View(Lambda, row, column, 1, 1));

            // Lambda_view = V_{j,\tau} = symmetrize( L^{-1}(\tau,x_k^{(j)}) . U_{j,k} )
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
      // if max_eigenvalue >= 0.
      const size_t num_eigvals(eigenvalues.Height());
      auto max_eigenvalue = eigenvalues(num_eigvals - 1, 0);
      ASSERT_EQUAL(max_eigenvalue, El::Max(eigenvalues),
                   "Eigenvalues were not sorted by El::HermitianEig()!");
      if(max_eigenvalue < 0)
        {
          PRINT_WARNING("block_", pvm_info.block_index,
                        ": x=", zero_values.at(zero_index),
                        ": negative max_eigenvalue=", max_eigenvalue,
                        " for Lambda matrix will be replaced with 0.");
          max_eigenvalue = 0;
        }
      {
        zeros.emplace_back(zero_values[zero_index]);
        auto &lambda(zeros.back().lambda);
        // lambdas = eigenvectors * sqrt(eigenvalues)
        // lambdas = v_{j,\tau} from Eq. (A.8)
        lambda = El::View(eigenvectors, 0, num_eigvals - 1, num_eigvals, 1);
        lambda *= El::Sqrt(max_eigenvalue);

        size_t row_column(0);
        for(size_t column(0); column != num_rows; ++column)
          for(size_t row(0); row <= column; ++row)
            {
              for(size_t index(0); index != matrix_block_size; ++index)
                {
                  error_matrix(index, row_column)
                    -= interpolation(index, zero_index) * lambda(row)
                       * lambda(column) * (row == column ? 1 : 2);
                }
              ++row_column;
            }

        // Set lambda = 1/sqrt(\chi) * v_{j,\tau}
        // With this definition, lambda does not change
        // if one adds reducedPrefactor != prefactor to PMP.json
        // NB: this is different from Python script and from (A.8) definition!
        lambda *= 1
                  / El::Sqrt(to_BigFloat(pvm_info.reduced_prefactor.evaluate(
                    to_Boost_Float(zeros.back().zero))));
      }
    }
  error = El::Sqrt(El::Dot(error_matrix, error_matrix));
}
