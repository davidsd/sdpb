#include <El.hpp>

El::BigFloat solve_LP(const El::Matrix<El::BigFloat> &A,
                      const El::Matrix<El::BigFloat> &b,
                      const El::Matrix<El::BigFloat> &c)
{
  El::Matrix<El::BigFloat> x, y, z;
  El::lp::direct::Ctrl<El::BigFloat> ctrl(false);
  // ctrl.mehrotraCtrl.print=true;
  // ctrl.mehrotraCtrl.maxIts=10;
  El::LP( A, b, c, x, y, z, ctrl );
  return (x(0)-x(1));
}
