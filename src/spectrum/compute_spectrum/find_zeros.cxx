#include "pmp/Polynomial_Vector_Matrix.hxx"
#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/assert.hxx"
#include "interpolate.hxx"

#include <El.hpp>

std::vector<El::BigFloat>
find_real_positive_minima_sorted(const Boost_Polynomial &polynomial);

namespace
{
  // For each block j, we want to build a matrix function that is equal to (c-B.y) at the sample points x_k
  // m_j_{r, s}(x_k) = (c - B.y) _{j, r, s, k}
  //
  // m_j_{r, s}(x)
  //   = p_{j, r, s}(x)*reduced_prefactor_j(x) * pv_j_r(x) * pv_j_s(x)
  //   p is polynomial, pv_j - preconditioning vector
  // How to build it:
  // - divide (c-B.y)_{j,r,s,k} by reduced_prefactor_j(x_k) * pv_j_r(x_k) * pv_j_s(x_k)
  // - for each {j,r,s}, build interpolating polynomial p_{j,r,s}(x), degree = (num_points - 1)
  Simple_Matrix<Boost_Polynomial> get_interpolated_polynomial_matrix(
    const El::Matrix<El::BigFloat> &c_minus_By_block,
    const Polynomial_Vector_Matrix &pvm)
  {
    const int height = pvm.polynomials.Height();
    const int width = pvm.polynomials.Width();
    ASSERT_EQUAL(height, width);
    const size_t num_points = pvm.sample_points.size();

    Simple_Matrix<Boost_Polynomial> interpolation_matrix(height, width);

    const bool has_preconditioning = pvm.preconditioning_vector.has_value();
    std::vector<std::vector<El::BigFloat>> pv_sampled(height);
    if(has_preconditioning)
      {
        for(int row = 0; row < height; ++row)
          {
            auto &samples = pv_sampled.at(row);
            samples.resize(num_points);
            for(size_t k = 0; k < num_points; ++k)
              {
                const auto &x_k = pvm.sample_points.at(k);
                samples.at(k)
                  = pvm.preconditioning_vector->at(row).evaluate(x_k);
              }
          }
      }

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
                  auto scale = pvm.reduced_sample_scalings.at(k);
                  if(has_preconditioning)
                    {
                      scale *= pv_sampled.at(i).at(k) * pv_sampled.at(j).at(k);
                    }
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

  El::BigFloat eval_determinant(
    const Simple_Matrix<Boost_Polynomial> &interpolated_poly_matrix,
    const Damped_Rational &reduced_prefactor,
    const std::optional<std::vector<Polynomial_Power_Product>>
      &preconditioning_vector,
    const Boost_Float &x)
  {
    const auto height = interpolated_poly_matrix.Height();
    const auto width = interpolated_poly_matrix.Width();
    ASSERT_EQUAL(height, width);

    std::optional<std::vector<Boost_Float>> pv;
    if(preconditioning_vector.has_value())
      {
        ASSERT_EQUAL(height, preconditioning_vector->size());
        pv->reserve(height);
        for(auto &ppp : preconditioning_vector.value())
          pv->emplace_back(ppp.evaluate(x));
      }
    const auto scale = reduced_prefactor.evaluate(x);

    El::Matrix<El::BigFloat> result(height, width);
    for(int i = 0; i < height; ++i)
      for(int j = 0; j < width; ++j)
        {
          Boost_Float value
            = interpolated_poly_matrix(i, j).evaluate(x) * scale;
          if(pv.has_value())
            value *= pv->at(i) * pv->at(j);
          result(i, j) = to_BigFloat(value);
        }
    return El::Determinant(result);
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
              const std::vector<El::BigFloat> &original_sample_points)
  {
    const int height = interpolated_poly_matrix.Height();
    const int width = interpolated_poly_matrix.Width();
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

    return interpolate(det_sample_points, det_samples);
  }

}

std::vector<El::BigFloat>
find_zeros(const El::Matrix<El::BigFloat> &c_minus_By_block,
           const Polynomial_Vector_Matrix &pvm, const El::BigFloat &threshold)
{
  ASSERT(threshold > 0, DEBUG_STRING(threshold));

  // Special case: constant constraint, isolated zero
  // In this case c-B.y is a constant matrix.
  // If it has a small eigenvalue, then we should report that we found a zero.
  // For convenience, we return zero at x=0
  // (in fact it doesn't depend on x, but spectrum.json format requires to specify zero value).
  // Note that all other zeros will be strictly positive.
  if(pvm.sample_points.size() == 1)
    {
      // El::HermitianEig modifies input matrix, so we have to make a copy
      auto block = c_minus_By_block;

      El::Matrix<El::BigFloat> eigenvalues;
      // Parameter tuning - copied from src/sdp_solve/SDP_Solver/run/step/step_length/min_eigenvalue.cxx
      El::HermitianEigCtrl<El::BigFloat> hermitian_eig_ctrl;
      hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.cutoff = block.Height() / 2 + 1;
      hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.secularCtrl.maxIterations
        = 16384;

      El::HermitianEig(El::UpperOrLowerNS::LOWER, block, eigenvalues,
                       hermitian_eig_ctrl);
      auto min_eigenvalue = El::Min(eigenvalues);
      ASSERT(min_eigenvalue > -threshold, "All eigenvalues must be positive!",
             DEBUG_STRING(min_eigenvalue), DEBUG_STRING(threshold));
      if(min_eigenvalue < threshold)
        return {0};
      return {};
    }

  // Regular case: polynomial matrix constraints.
  // First, we divide (c-B.y) by sample scalings, reduced_prefactor and preconditioning,
  // so that it corresponds to pure polynomials
  // and can be interpolated accordingly.
  // After that, we find minima of determinant.
  // Then we multiply by reduced_prefactor and preconditioning again,
  // and check if the minima are deep enough to be considered zeros.

  const auto interpolated_poly_matrix
    = get_interpolated_polynomial_matrix(c_minus_By_block, pvm);
  std::vector<El::BigFloat> minima = find_real_positive_minima_sorted(
    determinant(interpolated_poly_matrix, pvm.sample_points));

  if(minima.empty() || minima.back() > 0)
    {
      // We should always check x=0, even if MPSolve didn't find a root there
      minima.insert(minima.begin(), 0);
    }

  std::vector<El::BigFloat> zeros;

  // TODO compare min eigenvalues instead of determinants?
  const auto eval = [&](const El::BigFloat &x) {
    return eval_determinant(interpolated_poly_matrix, pvm.reduced_prefactor,
                            pvm.preconditioning_vector, to_Boost_Float(x));
  };

  for(size_t i = 0; i < minima.size(); ++i)
    {
      const auto &x = minima.at(i);
      const auto y = eval(x);

      bool is_zero = false;
      if(i == 0)
        {
          if(minima.size() > 1)
            {
              const auto x_right = get_midpoint(x, minima.at(i + 1));
              const auto y_right = eval(x_right);
              is_zero = y / y_right < threshold;
            }
          else
            {
              // This is a case of single minimum.
              // TODO: which points should we choose for comparison?
              // It's not obvious, choosing x/2 is not justified well.
              auto x_other = x / 2;
              if(x_other == El::BigFloat(0))
                {
                  // Special case: if x=0, we take the first nonzero sample point.
                  x_other = pvm.sample_points.at(0);
                  if(x_other == El::BigFloat(0))
                    x_other = pvm.sample_points.at(1);
                }
              ASSERT(x_other > 0);

              const auto y_other = eval(x_other);
              is_zero = y / y_other < threshold;
            }
        }
      else if(i + 1 == minima.size())
        {
          const auto x_left = get_midpoint(x, minima.at(i - 1));
          const auto y_left = eval(x_left);
          is_zero = y / y_left < threshold;
        }
      else
        {
          const auto x_left = get_midpoint(x, minima.at(i - 1));
          const auto y_left = eval(x_left);
          const auto x_right = get_midpoint(x, minima.at(i + 1));
          const auto y_right = eval(x_right);
          is_zero = y * y / y_left / y_right < threshold * threshold;
        }

      if(is_zero)
        zeros.emplace_back(x);
    }
  return zeros;
}
