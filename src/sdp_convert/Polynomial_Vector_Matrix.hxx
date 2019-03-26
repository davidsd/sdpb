#pragma once

#include "../Polynomial.hxx"

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
// is equivalent to a DualConstraintGroup
//
//   Tr(A_p Y) + (B y)_p = c_p
//
// A PolynomialVectorMatrix contains the data needed to construct this
// DualConstraintGroup:
//
class Polynomial_Vector_Matrix
{
public:
  int rows; // rows of M
  int cols; // cols of M

  // elements of M, in row-major order
  std::vector<std::vector<Polynomial>> elements;

  // A list of real numbers x_k (0 <= k <= degree(M)) at which to
  // sample M(x) to construct the v_{b,k}.
  std::vector<El::BigFloat> sample_points;

  // A list of real numbers s_k (0 <= k <= degree(M)) to scale M(x_k)
  // and the corresponding v_{b,k}.
  std::vector<El::BigFloat> sample_scalings;

  // bilinearBasis[m] = q_m(x) (0 <= m <= degree/2), where q_m is a
  // polynomial with degree deg(q_m) = m.
  std::vector<Polynomial> bilinear_basis;

  inline const std::vector<Polynomial> &elt(const int r, const int c) const
  {
    return elements[r + c * rows];
  }

  inline std::vector<Polynomial> &elt(const int r, const int c)
  {
    return elements[r + c * rows];
  }

  inline void clear()
  {
    elements.clear();
    sample_points.clear();
    sample_scalings.clear();
    bilinear_basis.clear();
  }
  
  // The maximal degree of any of the components P^{rs}_n(x).
  int64_t degree() const
  {
    int64_t d = 0;
    for(auto &e : elements)
      for(auto &p : e)
        {
          d = std::max(p.degree(), d);
        }
    return d;
  }
};
