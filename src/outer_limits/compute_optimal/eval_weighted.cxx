#include "poles_prefactor.hxx"
#include "power_prefactor.hxx"
#include "../../sdp_read.hxx"

El::BigFloat
eval_weighted(const Positive_Matrix_With_Prefactor &matrix,
              const El::BigFloat &x, const std::vector<El::BigFloat> &weights)
{
  const size_t matrix_dim(matrix.polynomials.size());
  El::Matrix<El::BigFloat> m(matrix_dim, matrix_dim);
  for(size_t row(0); row != matrix_dim; ++row)
    for(size_t column(0); column <= row; ++column)
      {
        auto &polys(matrix.polynomials.at(row).at(column));
        if(weights.size() != polys.size())
          {
            throw std::runtime_error("INTERNAL ERROR mismatch: "
                                     + std::to_string(weights.size()) + " "
                                     + std::to_string(polys.size()));
          }
        El::BigFloat element(0);
        for(size_t index(0); index != weights.size(); ++index)
          {
            element += weights[index] * polys[index](x);
          }
        m.Set(row, column, element);
      }

  // FIXME: Use the square of the matrix rather than the smallest
  // eigenvalue?  That would map to B^T B.
  
  El::Matrix<El::BigFloat> eigenvalues;
  /// There is a bug in El::HermitianEig when there is more than
  /// one level of recursion when computing eigenvalues.  One fix
  /// is to increase the cutoff so that there is no more than one
  /// level of recursion.

  /// An alternate workaround is to compute both eigenvalues and
  /// eigenvectors, but that seems to be significantly slower.
  El::HermitianEigCtrl<El::BigFloat> hermitian_eig_ctrl;
  hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.cutoff = matrix_dim / 2 + 1;

  /// The default number of iterations is 40.  That is sometimes
  /// not enough, so we bump it up significantly.
  hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.secularCtrl.maxIterations = 400;
  El::HermitianEig(El::UpperOrLowerNS::LOWER, m, eigenvalues,
                   hermitian_eig_ctrl);
  El::BigFloat result(El::Min(eigenvalues));

  if(!matrix.damped_rational.is_constant())
    {
      result *= poles_prefactor(matrix.damped_rational.poles, x)
                * power_prefactor(matrix.damped_rational.base, x);
    }
  return result;
}
