// The second step of the two-set encoding (pmp2sdp --encodingPoints), as a
// function of its own: the SDP of one block, written at its approximation
// points, in; the same block written at given encoding points out.
//
// Input: a Dual_Constraint_Group whose rows are, for each matrix entry (r,s)
// and approximation point x_k (k = 0..K-1), sigma_k R^{rs}_n(x_k), with
// sigma_k the reduced sample scalings (what pmp2sdp's default encoding
// writes: B and c), and encoding points y_j with scalings w_j (j = 0..K-1).
// Output, for each entry (r,s):
//
//   row j = w_j sum_k l_k(y_j) row_k / sigma_k,
//
// l_k the Lagrange polynomials of the x_k: w_j times the degree-(K-1)
// interpolant of R at the approximation points, evaluated at y_j, computed
// from the values at the x_k (barycentric weights), not from polynomial
// coefficients.  The bilinear bases are the polynomials orthonormal for
// sum_j w_j delta(x - y_j) (even) and sum_j w_j y_j delta(x - y_j) (odd),
// times sqrt(w_j) and sqrt(w_j y_j), made by the three-term recurrence
// (sample_bilinear_basis_by_recurrence).  The objective and normalization
// do not change.  The same call re-encodes an SDP that is already written
// (a refresh of the points at an iterate, a neighbour's points) given its
// rows; no blocks or polynomial coefficients are needed.
//
// Precision: the Lagrange step cancels.  The returned estimate is
//
//   max_j log10( max_n sum_k |l_k(y_j) V_kn| / max_n |sum_k l_k(y_j) V_kn| ),
//
// V_kn = row_k / sigma_k, the number of decimal digits of the input rows
// that the output rows lose relative to their largest entry (given inputs
// accurate to the working precision).  pmp2sdp runs this step at
// --encodingExtraPrecision bits above --precision.

#include "../Dual_Constraint_Group.hxx"
#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/assert.hxx"

#include <cmath>
#include <optional>

std::array<El::Matrix<El::BigFloat>, 2> sample_bilinear_basis_by_recurrence(
  const std::array<size_t, 2> &basis_sizes,
  const std::vector<El::BigFloat> &sample_points,
  const std::vector<El::BigFloat> &sample_scalings);

// The K x K matrix L(j, k) = l_k(y_j) of the Lagrange polynomials of the
// points x at the points y (barycentric form; exact 0/1 where y_j = x_k).
El::Matrix<El::BigFloat>
lagrange_matrix(const std::vector<El::BigFloat> &x,
                const std::vector<El::BigFloat> &y)
{
  const size_t K = x.size();
  ASSERT(K > 0);
  std::vector<El::BigFloat> weights(K);
  for(size_t k = 0; k < K; ++k)
    {
      El::BigFloat p = 1;
      for(size_t l = 0; l < K; ++l)
        if(l != k)
          {
            ASSERT(x[k] != x[l], "approximation points must be distinct",
                   DEBUG_STRING(k), DEBUG_STRING(l));
            p *= x[k] - x[l];
          }
      weights[k] = 1 / p;
    }
  El::Matrix<El::BigFloat> L(y.size(), K);
  for(size_t j = 0; j < y.size(); ++j)
    {
      std::optional<size_t> equal;
      El::BigFloat ell = 1;
      for(size_t k = 0; k < K; ++k)
        {
          if(y[j] == x[k])
            equal = k;
          ell *= y[j] - x[k];
        }
      for(size_t k = 0; k < K; ++k)
        {
          if(equal.has_value())
            L.Set(j, k, k == equal.value() ? 1 : 0);
          else
            L.Set(j, k, ell * weights[k] / (y[j] - x[k]));
        }
    }
  return L;
}

double reencode_dual_constraint_group(
  Dual_Constraint_Group &group,
  const std::vector<El::BigFloat> &approximation_points,
  const std::vector<El::BigFloat> &approximation_scalings,
  const std::vector<El::BigFloat> &encoding_points,
  const std::vector<El::BigFloat> &encoding_scalings,
  const std::array<size_t, 2> &basis_sizes)
{
  const size_t K = approximation_points.size();
  ASSERT_EQUAL(group.num_points, K);
  ASSERT_EQUAL(approximation_scalings.size(), K);
  ASSERT_EQUAL(encoding_points.size(), K,
               "the number of encoding points must be the number of "
               "approximation points");
  ASSERT_EQUAL(encoding_scalings.size(), K);
  const size_t num_entries = group.dim * (group.dim + 1) / 2;
  ASSERT_EQUAL(group.constraint_constants.size(), K * num_entries);
  ASSERT_EQUAL(group.constraint_matrix.Height(), K * num_entries);
  const size_t N = group.constraint_matrix.Width();

  const auto L = lagrange_matrix(approximation_points, encoding_points);
  El::Matrix<El::BigFloat> L_abs(K, K);
  for(size_t j = 0; j < K; ++j)
    for(size_t k = 0; k < K; ++k)
      L_abs.Set(j, k, El::Abs(L.CRef(j, k)));

  double digits_lost = 0;
  El::Matrix<El::BigFloat> V(K, N + 1), V_abs(K, N + 1), W(K, N + 1),
    W_abs(K, N + 1);
  const El::BigFloat zero(0);
  for(size_t e = 0; e < num_entries; ++e)
    {
      const size_t offset = e * K;
      for(size_t k = 0; k < K; ++k)
        {
          const El::BigFloat inv = 1 / approximation_scalings[k];
          const El::BigFloat c = group.constraint_constants[offset + k] * inv;
          V.Set(k, 0, c);
          V_abs.Set(k, 0, El::Abs(c));
          for(size_t n = 0; n < N; ++n)
            {
              const El::BigFloat b
                = group.constraint_matrix.CRef(offset + k, n) * inv;
              V.Set(k, n + 1, b);
              V_abs.Set(k, n + 1, El::Abs(b));
            }
        }
      El::Gemm(El::NORMAL, El::NORMAL, El::BigFloat(1), L, V, zero, W);
      El::Gemm(El::NORMAL, El::NORMAL, El::BigFloat(1), L_abs, V_abs, zero,
               W_abs);
      for(size_t j = 0; j < K; ++j)
        {
          El::BigFloat row_max = 0, row_abs_max = 0;
          for(size_t n = 0; n <= N; ++n)
            {
              row_max = std::max(row_max, El::Abs(W.CRef(j, n)));
              row_abs_max = std::max(row_abs_max, W_abs.CRef(j, n));
            }
          if(row_max > zero)
            digits_lost = std::max(
              digits_lost,
              static_cast<double>(log10(to_Boost_Float(row_abs_max / row_max))));
          const El::BigFloat &w = encoding_scalings[j];
          group.constraint_constants[offset + j] = w * W.CRef(j, 0);
          for(size_t n = 0; n < N; ++n)
            group.constraint_matrix.Set(offset + j, n, w * W.CRef(j, n + 1));
        }
    }

  group.bilinear_bases = sample_bilinear_basis_by_recurrence(
    basis_sizes, encoding_points, encoding_scalings);
  return digits_lost;
}
