#include "outer_limits/Function.hxx"
#include <El.hpp>
#include <vector>

El::BigFloat
eval_summed(const El::BigFloat &epsilon, const El::BigFloat &infinity,
            const std::vector<std::vector<Function>> &summed_functions,
            const El::BigFloat &x)
{
  const size_t matrix_dim(summed_functions.size());
  El::Matrix<El::BigFloat> m(matrix_dim, matrix_dim);
  for(size_t row(0); row != matrix_dim; ++row)
    for(size_t column(0); column <= row; ++column)
      {
        m.Set(row, column,
              summed_functions.at(row).at(column).eval(epsilon, infinity, x));
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
  return El::Min(eigenvalues);
}
