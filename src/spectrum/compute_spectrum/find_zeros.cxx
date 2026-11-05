#include "interpolate.hxx"
#include "pmp/PMP_Info.hxx"
#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <El.hpp>

std::vector<El::BigFloat>
find_real_positive_minima_sorted(const Boost_Polynomial &polynomial,
                                 const El::BigFloat &min_zero_distance,
                                 Timers &timers);

namespace
{
  // For each block j, we want to build a matrix function that is equal to (c-B.y) at the sample points x_k
  // m_j_{r, s}(x_k) = (c - B.y) _{j, r, s, k}
  //
  // m_j_{r, s}(x)
  //   = p_{j, r, s}(x)*reduced_prefactor_j(x)
  //   p is polynomial
  // How to build it:
  // - divide (c-B.y)_{j,r,s,k} by reduced_prefactor_j(x_k)
  // - for each {j,r,s}, build interpolating polynomial p_{j,r,s}(x), degree = (num_points - 1)
  Simple_Matrix<Boost_Polynomial> get_interpolated_polynomial_matrix(
    const El::Matrix<El::BigFloat> &c_minus_By_block, const PVM_Info &pvm,
    Timers &timers)
  {
    Scoped_Timer timer(timers, "interpolate");
    const auto height = pvm.dim;
    const auto width = height;
    const size_t num_points = pvm.sample_points.size();

    Simple_Matrix<Boost_Polynomial> interpolation_matrix(height, width);

    const auto lagrange_basis = get_lagrange_basis(pvm.sample_points);

    ASSERT_EQUAL(c_minus_By_block.Height(),
                 height * (height + 1) / 2 * num_points);
    ASSERT_EQUAL(c_minus_By_block.Width(), 1);
    {
      int rsk_index = 0;
      for(int i = 0; i < height; ++i)
        {
          for(int j = 0; j <= i; ++j)
            {
              std::vector<El::BigFloat> ys;
              for(size_t k = 0; k < num_points; ++k)
                {
                  const auto &scale = pvm.reduced_sample_scalings.at(k);
                  ys.push_back(c_minus_By_block.CRef(rsk_index) / scale);

                  ++rsk_index;
                }

              interpolation_matrix(i, j) = interpolate(lagrange_basis, ys);
              // Symmetrize
              if(i != j)
                interpolation_matrix(j, i) = interpolation_matrix(i, j);
            }
        }
    }
    return interpolation_matrix;
  }

  Boost_Float eval_determinant(
    const Simple_Matrix<Boost_Polynomial> &interpolated_poly_matrix,
    const Damped_Rational &reduced_prefactor, const Boost_Float &x)
  {
    const auto height = interpolated_poly_matrix.Height();
    const auto width = interpolated_poly_matrix.Width();
    ASSERT_EQUAL(height, width);

    El::Matrix<El::BigFloat> poly_values(height, width);
    for(int i = 0; i < height; ++i)
      for(int j = 0; j < width; ++j)
        {
          poly_values(i, j)
            = to_BigFloat(interpolated_poly_matrix(i, j).evaluate(x));
        }

    const auto det = to_Boost_Float(El::Determinant(poly_values));
    // Multiply each element by chi, so that det is multiplied by chi^N
    const auto chi = reduced_prefactor.evaluate(x);
    return det * pow(chi, height);
  }

  template <class T> T get_midpoint(const T &a, const T &b)
  {
    ASSERT(a != b, DEBUG_STRING(a), "Points should be different!");

    // Harmonic mean will return 0, so let's use arithmetic mean.
    if(a == T(0) || b == T(0))
      return (a + b) / 2;

    // Rajeev argued that harmonic mean works better in his tests.
    return 2 * a * b / (a + b);
  }

  Boost_Polynomial
  determinant(const Simple_Matrix<Boost_Polynomial> &interpolated_poly_matrix,
              const std::vector<El::BigFloat> &original_sample_points,
              Timers &timers)
  {
    Scoped_Timer timer(timers, "determinant");
    const auto height = interpolated_poly_matrix.Height();
    const auto width = interpolated_poly_matrix.Width();
    ASSERT_EQUAL(height, width);

    // Trivial case: determinant of a 1x1 matrix
    if(height == 1)
      return interpolated_poly_matrix(0, 0);

    ASSERT(!original_sample_points.empty());
    // Max degree of element
    const int element_degree = original_sample_points.size() - 1;
    // Degree of determinant
    const int det_degree = element_degree * height;

    // Number of sampling points to define polynomial determinant
    const int det_num_points = det_degree + 1;

    // Sample points for resulting determinant polynomial
    std::vector<El::BigFloat> det_sample_points;
    det_sample_points.reserve(det_num_points);

    // Add (height - 1) evenly spaced points between each pair of adjacent sample points
    // In total, the number of points is
    // num_points + (height - 1) * (num_points - 1) = height * degree + 1 = det_degree + 1
    for(size_t i = 0; i + 1 < original_sample_points.size(); i++)
      {
        const auto &x = original_sample_points.at(i);
        const auto &x_next = original_sample_points.at(i + 1);
        const El::BigFloat delta = (x_next - x) / height;
        for(size_t k = 0; k < height; ++k)
          {
            det_sample_points.emplace_back(x + delta * k);
          }
      }
    det_sample_points.emplace_back(original_sample_points.back());

    ASSERT_EQUAL(det_sample_points.size(), det_num_points);

    std::vector<El::BigFloat> det_samples;
    det_samples.reserve(det_num_points);

    {
      El::Matrix<El::BigFloat> m(height, width);
      for(const auto &x : det_sample_points)
        {
          const Boost_Float xx = to_Boost_Float(x);
          for(int i = 0; i < height; ++i)
            for(int j = 0; j < width; ++j)
              {
                m(i, j)
                  = to_BigFloat(interpolated_poly_matrix(i, j).evaluate(xx));
              }

          det_samples.emplace_back(El::Determinant(m));
        }
    }

    Scoped_Timer interpolate_timer(timers, "interpolate");
    return interpolate(det_sample_points, det_samples);
  }

}

std::vector<El::BigFloat>
find_zeros(const El::Matrix<El::BigFloat> &c_minus_By_block,
           const PVM_Info &pvm, const Boost_Float &threshold,
           const El::BigFloat &max_zero, const El::BigFloat &min_zero_distance,
           Timers &timers)
{
  Scoped_Timer timer(timers, "find_zeros");
  ASSERT(threshold > 0, DEBUG_STRING(threshold));

  // Special case: constant constraint, isolated zero
  // In this case c-B.y is a constant matrix.
  // If it has a small eigenvalue, then we should report that we found a zero.
  // For convenience, we return zero at x=0
  // (in fact it doesn't depend on x, but spectrum.json format requires to specify zero value).
  // Note that all other zeros will be strictly positive.
  if(pvm.sample_points.size() == 1)
    {
      Scoped_Timer const_timer(timers, "constant_constraint");

      // c - B.y is a vector, we need to reshape it into a (dim x dim) matrix.
      // We reuse get_interpolated_polynomial_matrix() to build it.
      // In this case the interpolation matrix contains 0-degree polynomials:
      // interpolated_poly_matrix[r,s] = (c - B.y)[rs,1] / scale
      // where scale is reduced_sample_scaling at x_0.
      const auto interpolated_poly_matrix
        = get_interpolated_polynomial_matrix(c_minus_By_block, pvm, timers);
      const auto dim = pvm.dim;
      ASSERT_EQUAL(interpolated_poly_matrix.Height(), dim);
      ASSERT_EQUAL(interpolated_poly_matrix.Width(), dim);

      El::Matrix<El::BigFloat> block(dim, dim);
      const auto &x = pvm.sample_points.front();
      ASSERT_EQUAL(pvm.reduced_sample_scalings.size(), 1);
      const auto &scale = pvm.reduced_sample_scalings.front();
      for(int i = 0; i < dim; ++i)
        {
          for(int j = 0; j < dim; ++j)
            {
              const auto value
                = interpolated_poly_matrix(i, j).evaluate(to_Boost_Float(x));
              // NB: restore scaling removed in get_interpolated_polynomial_matrix()
              block.Set(i, j, to_BigFloat(value) * scale);
            }
        }

      El::Matrix<El::BigFloat> eigenvalues;
      // Parameter tuning - copied from src/sdp_solve/SDP_Solver/run/step/step_length/min_eigenvalue.cxx
      El::HermitianEigCtrl<El::BigFloat> hermitian_eig_ctrl;
      hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.cutoff = block.Height() / 2 + 1;
      hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.secularCtrl.maxIterations
        = 16384;

      El::HermitianEig(El::UpperOrLowerNS::LOWER, block, eigenvalues,
                       hermitian_eig_ctrl);
      auto min_eigenvalue = to_Boost_Float(El::Min(eigenvalues));
      ASSERT(min_eigenvalue > -threshold, "All eigenvalues must be positive!",
             DEBUG_STRING(min_eigenvalue), DEBUG_STRING(threshold));
      if(min_eigenvalue < threshold)
        return {0};
      return {};
    }

  // Regular case: polynomial matrix constraints.
  // First, we divide (c-B.y) by sample scalings and reduced_prefactor,
  // so that it corresponds to pure polynomials
  // and can be interpolated accordingly.
  // After that, we find minima of determinant.
  // Then we multiply by reduced_prefactor,
  // and check if the minima are deep enough to be considered zeros.

  const auto interpolated_poly_matrix
    = get_interpolated_polynomial_matrix(c_minus_By_block, pvm, timers);
  const auto det
    = determinant(interpolated_poly_matrix, pvm.sample_points, timers);
  std::vector<El::BigFloat> minima;
  for(auto &x :
      find_real_positive_minima_sorted(det, min_zero_distance, timers))
    {
      // Remove large zeros
      if(max_zero > 0 && x > max_zero)
        {
          PRINT_WARNING("block_", pvm.block_index,
                        ": ignore large zero at x=", x);
          break;
        }
      minima.push_back(x);
    }

  if(minima.empty() || minima.front() > 0)
    {
      // We should always check x=0, even if MPSolve didn't find a root there
      minima.insert(minima.begin(), 0);
    }

  std::vector<El::BigFloat> zeros;

  // TODO compare min eigenvalues instead of determinants?
  const auto eval = [&](const El::BigFloat &x) {
    return eval_determinant(interpolated_poly_matrix, pvm.reduced_prefactor,
                            to_Boost_Float(x));
  };

  const auto is_zero_one_sided
    = [&](const El::BigFloat &x, const El::BigFloat &x_other) {
        const auto y = eval(x);
        const auto y_other = eval(x_other);
        const auto ratio = y / y_other;
        ASSERT(!isnan(ratio), "Cannot check for zero: det(c-By)=", y,
               " at x=", x, ", det(c-By)=", y_other, " at x=", x_other);
        return ratio < threshold;
      };

  // TODO should it simply return
  // is_zero_one_sided(x, x_left) && is_zero_one_sided(x, x_right)?
  const auto is_zero_two_sided
    = [&](const El::BigFloat &x, const El::BigFloat &x_left,
          const El::BigFloat &x_right) {
        const auto y = eval(x);
        const auto y_left = eval(x_left);
        const auto y_right = eval(x_right);
        const auto ratio_squared = y * y / y_left / y_right;
        ASSERT(!isnan(ratio_squared), "Cannot check for zero: det(c-By)=", y,
               " at x=", x, ", det(c-By)=", y_left, " at x=", x_left,
               ", det(c-By)=", y_right, " at x=", x_right);
        return ratio_squared < threshold * threshold;
      };

  Scoped_Timer check_minima_timer(timers, "check_minima");
  for(size_t i = 0; i < minima.size(); ++i)
    {
      const auto &x = minima.at(i);
      // TODO: what if y (or its neighbor) is infinite or NaN?
      // This is possible e.g. if x=0 and prefactor has pole at x=0.

      const bool is_zero = [&] {
        // First zero candidate
        if(i == 0)
          {
            if(minima.size() > 1)
              {
                const auto x_right = get_midpoint(x, minima.at(i + 1));
                return is_zero_one_sided(x, x_right);
              }
            // This is a case of single minimum.
            // TODO: which points should we choose for comparison?
            // It's not obvious, choosing x/2 is not justified well.
            auto x_other = x / 2;
            if(x_other == El::BigFloat(0))
              {
                // Special case: if x=0, we take the first nonzero sample point.
                // Note that we have at least two sample points:
                // constant constraints are handled separately (see above).
                x_other = pvm.sample_points.at(0);
                if(x_other == El::BigFloat(0))
                  x_other = pvm.sample_points.at(1);
              }
            ASSERT(x_other > 0);
            return is_zero_one_sided(x, x_other);
          }
        // Last zero candidate
        if(i + 1 == minima.size())
          {
            const auto x_left = get_midpoint(x, minima.at(i - 1));
            return is_zero_one_sided(x, x_left);
          }
        // Regular case, check left and right neighbors
        const auto x_left = get_midpoint(x, minima.at(i - 1));
        const auto x_right = get_midpoint(x, minima.at(i + 1));
        return is_zero_two_sided(x, x_left, x_right);
      }();

      if(is_zero)
        zeros.emplace_back(x);
    }
  return zeros;
}
