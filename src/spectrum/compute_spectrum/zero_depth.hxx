#pragma once

#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/assert.hxx"

#include <El.hpp>
#include <limits>
#include <vector>

// Local zero test.
//
// Let M(x) = (c-B.y)(x) / chi(x) be the polynomial matrix of a block
// (prefactor chi removed), and P(x) = det M(x),
// so that det(c-B.y) = P(x) * chi(x)^dim with a positive prefactor chi(x).
// A candidate x0 is either a local minimum of P or the boundary point x0 = 0.
//
// A zero of the functional is a dip of P that is narrow compared to the
// distance L to the nearest other structure of P. We measure it with the
// logarithmic derivatives of P at x0, computed exactly from M by Jacobi's
// formula (no determinant polynomial, no finite differences):
//
//   P'/P  = tr(M^-1 M'),
//   P''/P = tr(M^-1 M'') - tr((M^-1 M')^2) + (P'/P)^2,
//
//   minimum:  depth = 1 / max(L |P'/P|, L^2 (P''/P) / 2)   (infinite if P'' <= 0)
//   boundary: depth = 1 / (L P'/P)                          (infinite if P' <= 0)
//
// For a quadratic dip P ~ c ((x-x0)^2 + w^2) this is depth = (w/L)^2:
// the zero is located to relative accuracy sqrt(depth) of L. If several
// eigenvalues of M dip at x0, the depth is that of the narrowest dip
// (not their product). At x0 = 0 a zero is a simple zero just left of 0,
// at distance delta: depth = delta/L. The L |P'/P| term rejects candidates
// that are not actually minima of P.
//
// L is the distance from x0 to the nearest other critical point of P in
// (0, inf); for a minimum with no critical point on its left, the boundary
// x = 0 also counts. If P has no other critical point, L = fallback_length.
// L is capped at fallback_length (the largest sample point): MPSolve can
// return a spurious huge root of P' (1e9..1e100) when P is monotone, which
// would otherwise make every boundary candidate look deep.
//
// The prefactor does not enter: it is smooth and positive, so it changes
// neither the location nor the width of a dip, while its exponential decay
// is what makes comparisons with distant points unreliable.
//
// Known misses:
// - a double zero at x0 > 0 with no critical point on its left gets L = x0,
//   so depth = (w/x0)^2: it is missed when x0 < w/sqrt(threshold), i.e. a
//   zero closer to x = 0 than ~1/sqrt(threshold) of its own widths (the
//   boundary candidate then has P'(0) < 0, depth infinite);
// - a flat (e.g. quartic) dip has P''(x0) ~ 0, depth infinite.
//
// The candidate is a zero iff depth < threshold.
// P(x0) <= 0 gives depth <= 0, i.e. a zero.

inline El::BigFloat
zero_local_length(const El::BigFloat &x, const bool is_boundary,
                  const std::vector<El::BigFloat> &sorted_critical_points,
                  const El::BigFloat &min_distance,
                  const El::BigFloat &fallback_length)
{
  El::BigFloat length = -1;
  bool has_left = false;
  for(const auto &c : sorted_critical_points)
    {
      if(c <= 0)
        continue;
      const El::BigFloat distance = El::Abs(c - x);
      if(distance <= min_distance)
        continue;
      if(c < x)
        has_left = true;
      if(length < 0 || distance < length)
        length = distance;
    }
  if(!is_boundary && !has_left && x > 0 && (length < 0 || x < length))
    length = x;
  if(length < 0 || length > fallback_length)
    length = fallback_length;
  ASSERT(length > 0, DEBUG_STRING(x), DEBUG_STRING(length));
  return length;
}

// M, dM, d2M: values of M(x), M'(x), M''(x) at the candidate.
inline Boost_Float zero_local_depth(const El::Matrix<El::BigFloat> &M,
                                    const El::Matrix<El::BigFloat> &dM,
                                    const El::Matrix<El::BigFloat> &d2M,
                                    const El::BigFloat &length,
                                    const bool is_boundary)
{
  const int dim = M.Height();
  ASSERT(M.Width() == dim && dM.Height() == dim && dM.Width() == dim
         && d2M.Height() == dim && d2M.Width() == dim);
  ASSERT(length > 0);
  const Boost_Float infinity = std::numeric_limits<Boost_Float>::infinity();

  const El::BigFloat det = El::Determinant(M);
  if(det <= El::BigFloat(0))
    return to_Boost_Float(det); // zero (or below zero)

  El::Matrix<El::BigFloat> M_inv(M);
  El::Inverse(M_inv);
  El::Matrix<El::BigFloat> A, B, AA;
  El::Gemm(El::NORMAL, El::NORMAL, El::BigFloat(1), M_inv, dM, A);
  El::Gemm(El::NORMAL, El::NORMAL, El::BigFloat(1), M_inv, d2M, B);
  El::Gemm(El::NORMAL, El::NORMAL, El::BigFloat(1), A, A, AA);
  El::BigFloat log_d1 = 0, log_d2 = 0;
  for(int i = 0; i < dim; ++i)
    {
      log_d1 += A(i, i);
      log_d2 += B(i, i) - AA(i, i);
    }
  // P'/P and P''/P
  const Boost_Float d1 = to_Boost_Float(log_d1);
  const Boost_Float d2 = to_Boost_Float(log_d2 + log_d1 * log_d1);
  const Boost_Float L = to_Boost_Float(length);

  if(is_boundary)
    return d1 > 0 ? Boost_Float(1 / (L * d1)) : infinity;
  if(d2 <= 0)
    return infinity;
  return 1 / std::max(Boost_Float(L * abs(d1)), Boost_Float(L * L * d2 / 2));
}
