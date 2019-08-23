// Given a vector of polynomials {q_0(x), q_1(x), ..., q_n(x)} of
// degree deg q_m(x) = m, a list of numSamples points x_k and scaling
// factors s_k, form the (maxDegree+1) x numSamples Matrix
//
//   {{ \sqrt(s_0) q_0(x_0), ..., \sqrt(s_K) q_0(x_K) },
//    { \sqrt(s_0) q_1(x_0), ..., \sqrt(s_K) q_1(x_K) },
//    ...
//    { \sqrt(s_0) q_M(x_0), ..., \sqrt(s_K) q_M(x_K) }}
//
// where maxDegree = M and numSamples = K+1.
//
// Input:
// - maxDegree: the maximal degree of q_m(x) to include
// - numSamples: number of sample points x_k
// - bilinearBasis: the vector {q_0(x), q_1(x), ..., q_n(x)}
// - samplePoints: the points {x_0, x_1, ... }
// - sampleScalings: the scale factors {s_0, s_1, ... }
//


#include "../../../sdp2input/write_output/bilinear_basis/factorial.hxx"

#include "../../../sdp2input/Boost_Float.hxx"
#include <boost/math/special_functions.hpp>

template Boost_Float
boost::math::chebyshev_clenshaw_recurrence<Boost_Float, Boost_Float>(
  const Boost_Float *c, size_t length, const Boost_Float &x);

El::Matrix<El::BigFloat>
sample_bilinear_basis(const int max_degree, const int num_samples,
                      const std::vector<Boost_Float> &sample_points,
                      const std::vector<Boost_Float> &sample_scalings)
{
  Boost_Float::default_precision(El::gmp::Precision() * log(2) / log(10));

  El::Matrix<El::BigFloat> b(max_degree + 1, num_samples);
  for(int k = 0; k < num_samples; k++)
    {
      Boost_Float x(sample_points[k]);
      Boost_Float scale(sqrt(sample_scalings[k]));

      for(int i = 0; i <= max_degree; i++)
        {
          std::vector<Boost_Float> cheb_polynomial(i + 1, Boost_Float(0));
          cheb_polynomial.back() = 1;
          Boost_Float b_boost(
            scale
            * boost::math::chebyshev_clenshaw_recurrence(
                cheb_polynomial.data(), cheb_polynomial.size(), x));
          b.Set(i, k, El::BigFloat(to_string(b_boost)));
        }
    }
  return b;
}
