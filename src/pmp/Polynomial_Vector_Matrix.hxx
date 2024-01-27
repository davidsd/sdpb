#pragma once

#include "Polynomial.hxx"
#include "Damped_Rational.hxx"

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
  El::Matrix<Polynomial_Vector> polynomials;

  // A list of real numbers x_k (0 <= k <= degree(M)) at which to
  // sample M(x) to construct the v_{b,k}.
  std::vector<El::BigFloat> sample_points;

  // A list of real numbers s_k (0 <= k <= degree(M)) to scale M(x_k)
  // and the corresponding v_{b,k}.
  std::vector<El::BigFloat> sample_scalings;

  // bilinearBasis[m] = q_m(x) (0 <= m <= degree/2), where q_m is a
  // polynomial with degree deg(q_m) = m.
  Polynomial_Vector bilinear_basis;

  Polynomial_Vector_Matrix(
    const El::Matrix<Polynomial_Vector> &polynomials,
    const std::optional<Damped_Rational> &prefactor_opt,
    const std::optional<std::vector<El::BigFloat>> &sample_points_opt,
    const std::optional<std::vector<El::BigFloat>> &sample_scalings_opt,
    const std::optional<Polynomial_Vector>
      &bilinear_basis_opt) noexcept(false);

  void validate(int64_t max_degree) const;
};
