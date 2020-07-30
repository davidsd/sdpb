#include <El.hpp>

void solve_LP(const El::Matrix<El::BigFloat> &A,
              const El::Matrix<El::BigFloat> &b,
              const El::Matrix<El::BigFloat> &c,
              std::vector<El::BigFloat> &weights)
{
  El::Matrix<El::BigFloat> x, y, z;
  El::lp::direct::Ctrl<El::BigFloat> ctrl(false);
  El::LP(A, b, c, x, y, z, ctrl);

  // W = w_+ - w_-
  for(size_t weight_index(0); weight_index != weights.size(); ++weight_index)
    {
      weights[weight_index]
        = x(2 * weight_index, 0) - x(2 * weight_index + 1, 0);
    }
}
