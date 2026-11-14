// --stepLengthDouble: see step_length_double.hxx
#include "step_length_double.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace
{
  // x = d 2^ex with 0.5 <= |d| < 1
  El::Int exponent_of(const El::BigFloat &x)
  {
    long ex;
    mpf_get_d_2exp(&ex, x.gmp_float.get_mpf_t());
    return ex;
  }

  // x 2^{-e} as a double (0 below double's range)
  double scaled_double(const El::BigFloat &x, const El::Int e)
  {
    long ex;
    const double d = mpf_get_d_2exp(&ex, x.gmp_float.get_mpf_t());
    if(d == 0)
      return 0;
    const El::Int s = ex - e;
    if(s < -1100)
      return 0;
    return std::ldexp(d, static_cast<int>(s));
  }
}

Double_Min_Eig double_min_eigenvalue(const El::DistMatrix<El::BigFloat> &block,
                                     const double gamma, const double tol)
{
  Double_Min_Eig out;
  const El::Int n = block.Height();
  if(n == 0)
    return out;
  const El::Int none = std::numeric_limits<El::Int>::min();

  // the block's scale, over its grid
  const El::Matrix<El::BigFloat> &local = block.LockedMatrix();
  El::Int e = none;
  for(El::Int j = 0; j < local.Width(); ++j)
    for(El::Int i = 0; i < local.Height(); ++i)
      if(local.CRef(i, j).Sign() != 0)
        e = std::max(e, exponent_of(local.CRef(i, j)));
  if(block.Grid().Size() > 1)
    e = El::mpi::AllReduce(e, El::mpi::MAX, block.DistComm());
  if(e == none)
    return out; // zero block

  // A 2^{-e} in double, gathered on every rank of the grid
  El::DistMatrix<double> A_d(block.Grid());
  A_d.AlignWith(block);
  A_d.Resize(n, n);
  El::Matrix<double> &A_d_local = A_d.Matrix();
  for(El::Int j = 0; j < local.Width(); ++j)
    for(El::Int i = 0; i < local.Height(); ++i)
      A_d_local(i, j) = scaled_double(local.CRef(i, j), e);
  const El::DistMatrix<double, El::STAR, El::STAR> A_ss(A_d);
  const El::Matrix<double> &p = A_ss.LockedMatrix();

  std::vector<double> A(size_t(n) * n), w(n);
  double fro = 0;
  for(El::Int j = 0; j < n; ++j)
    for(El::Int i = 0; i < n; ++i)
      {
        const double v = 0.5 * (p.CRef(i, j) + p.CRef(j, i));
        A[size_t(j) * n + i] = v;
        fro += v * v;
      }
  El::lapack::HermitianEig('L', n, A.data(), n, w.data());
  // ascending; take the minimum anyway
  const double lambda = *std::min_element(w.begin(), w.end());
  const double EPS = std::ldexp(1.0, -52);
  const double bound = (n + 2) * EPS * std::sqrt(fro);

  using ld = long double;
  auto step_from = [&](const ld lam) -> ld {
    return lam > -ld(gamma) ? ld(1) : -ld(gamma) / lam;
  };
  const ld l = std::ldexp(ld(lambda), int(e));
  const ld err = std::ldexp(ld(bound), int(e));
  if(step_from(l - err) >= (1 - ld(tol)) * step_from(l + err))
    {
      out.status = Double_Min_Eig_Status::accepted;
      out.value = El::BigFloat(lambda);
      if(e > 0)
        mpf_mul_2exp(out.value.gmp_float.get_mpf_t(),
                     out.value.gmp_float.get_mpf_t(), e);
      else if(e < 0)
        mpf_div_2exp(out.value.gmp_float.get_mpf_t(),
                     out.value.gmp_float.get_mpf_t(), -e);
    }
  else
    out.status = Double_Min_Eig_Status::flagged;
  return out;
}
