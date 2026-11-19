#pragma once

#include "Polynomial.hxx"
#include "sdpb_util/Damped_Rational.hxx"
#include "sdpb_util/Simple_Matrix.hxx"

#include <optional>

// Let M(x) be a matrix whose entries are vectors of polynomials:
//
//   M(x) = ( \vec P^{00}(x) ... \vec P^{m0}(x) )
//          ( ...                               )
//          ( \vec P^{0n}(x) ... \vec P^{mn}(x) )
//
// where each vector has length N+1:
//
//   \vec P^{rs}(x) = (P^{rs}_{-1}(x), P^{rs}_0, ... , P^{rs}_{N-1}(x))
//
// Consider a vector y = (y_0, ..., y_{N-1}) of length N, and let
// (1,y) denote the vector of length N+1 whose components are 1,
// followed by the components of y.  As explained in the manual, the
// constraint
//
//   (1,y) . M(x) is positive semidefinite
//
// is equivalent to a Dual_Constraint_Group
//
//   Tr(A_p Y) + (B y)_p = c_p
//
// A Polynomial_Vector_Matrix contains the data needed to construct this
// Dual_Constraint_Group.
// It also may contain Damped_Rational prefactor
struct Polynomial_Vector_Matrix
{
  // TODO: This is actually a symmetric matrix (dim x dim) of a vector
  // of polynomials.  So we should use a proper symmetric matrix class.
  Simple_Matrix<Polynomial_Vector> polynomials;

  Damped_Rational prefactor;
  Damped_Rational reduced_prefactor;

  // A list of real numbers x_k (0 <= k <= degree(M)) at which to
  // sample M(x) to construct the v_{b,k}.
  std::vector<El::BigFloat> sample_points;

  // A list of real numbers s_k (0 <= k <= degree(M)) to scale M(x_k)
  // and the corresponding v_{b,k}.
  std::vector<El::BigFloat> sample_scalings;
  std::vector<El::BigFloat> reduced_sample_scalings;

  // bilinearBasis[m] = q_m(x) (0 <= m <= degree/2), where q_m is a
  // polynomial with degree deg(q_m) = m.
  std::array<Polynomial_Vector, 2> bilinear_basis;

  // Two-set encoding (pmp2sdp --encodingPoints, see
  // convert/encoding_points.cxx).  If not empty, the block's polynomials are
  // still defined at sample_points (the approximation points), and the SDP
  // is written at these encoding points with these scalings (for B, c and
  // the bilinear bases alike), the bases made by recurrence.  Empty
  // (default): the SDP is written at the sample points.
  std::vector<El::BigFloat> encoding_points;
  std::vector<El::BigFloat> encoding_scalings;

  Polynomial_Vector_Matrix(
    const Simple_Matrix<Polynomial_Vector> &polynomials,
    const std::optional<Damped_Rational> &prefactor_opt,
    const std::optional<Damped_Rational> &reduced_prefactor_opt,
    const std::optional<int64_t> &max_num_poles_opt,
    const std::optional<std::vector<El::BigFloat>> &sample_points_opt,
    const std::optional<std::vector<El::BigFloat>> &sample_scalings_opt,
    const std::optional<std::vector<El::BigFloat>> &reduced_sample_scalings_opt,
    const std::optional<std::array<Polynomial_Vector, 2>>
      &bilinear_basis_opt) noexcept(false);

  void validate(int64_t num_points) const;
};

// Two-set encoding (pmp2sdp --encodingPoints; convert/encoding_points.cxx).
// The scale a = scale_factor * median |kept pole| (10 * scale_factor without
// poles).
El::BigFloat encoding_points_scale(const Damped_Rational &reduced_prefactor,
                                   const El::BigFloat &scale_factor);
// The functionalMatched rule: the points and weights of the quadrature that
// is exact to degree num_points - 1 for dx / (sqrt(x) g(x)), with g made of
// the kept poles and the scale a.
void functional_matched_encoding_points(
  const Damped_Rational &reduced_prefactor, size_t num_points,
  const El::BigFloat &scale, std::vector<El::BigFloat> &points,
  std::vector<El::BigFloat> &scalings);
// The chebyshev rule: Chebyshev points of x / (x + a), with the scalings
// 1 / (prod_i (x - p_i) (x + a)^(num_points - 1 - #poles)).
void chebyshev_encoding_points(const Damped_Rational &reduced_prefactor,
                               size_t num_points, const El::BigFloat &scale,
                               std::vector<El::BigFloat> &points,
                               std::vector<El::BigFloat> &scalings);
// Check and set m.encoding_points and m.encoding_scalings.
void set_encoding_points(Polynomial_Vector_Matrix &m,
                         std::vector<El::BigFloat> points,
                         std::vector<El::BigFloat> scalings);
