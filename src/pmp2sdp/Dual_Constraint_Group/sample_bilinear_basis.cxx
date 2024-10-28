// Given a vector of polynomials {q_0(x), q_1(x), ..., q_n(x)}
// of degree deg q_m(x) = m,
// a list of points {x_0 .. x_K}
// and scaling factors {s_0 .. s_K}, form the Matrix
//
//   {{ \sqrt(s_0) q_0(x_0), ..., \sqrt(s_K) q_0(x_K) },
//    { \sqrt(s_0) q_1(x_0), ..., \sqrt(s_K) q_1(x_K) },
//    ...
//    { \sqrt(s_0) q_n(x_0), ..., \sqrt(s_K) q_n(x_K) }}
//
// Input:
// - bilinearBasis: the vector {q_0(x), q_1(x), ..., q_n(x)}
// - samplePoints: the points {x_0, x_1, ... x_K}
// - sampleScalings: the scale factors {s_0, s_1, ... s_K}

#include "pmp/Polynomial.hxx"

El::Matrix<El::BigFloat>
sample_bilinear_basis(const Polynomial_Vector &bilinearBasis,
                      const std::vector<El::BigFloat> &samplePoints,
                      const std::vector<El::BigFloat> &sampleScalings)
{
  El::Matrix<El::BigFloat> b(bilinearBasis.size(), samplePoints.size());
  for(int k = 0; k < samplePoints.size(); k++)
    {
      const auto &x = samplePoints.at(k);
      const auto scale = Sqrt(sampleScalings.at(k));
      for(int i = 0; i < bilinearBasis.size(); i++)
        {
          b.Set(i, k, scale * bilinearBasis.at(i)(x));
        }
    }
  return b;
}

std::array<El::Matrix<El::BigFloat>, 2>
sample_bilinear_basis(const std::array<Polynomial_Vector, 2> &bilinear_basis,
                      const std::vector<El::BigFloat> &sample_points,
                      const std::vector<El::BigFloat> &sample_scalings)
{
  // The matrix Y has two blocks Y_1, Y_2.  The bilinear_bases for the
  // constraint matrices A_p are given by sampling the following
  // vectors for each block:
  //
  //   Y_1: {q_0(x), ..., q_delta1(x)}
  //   Y_2: {\sqrt(x) q_0(x), ..., \sqrt(x) q_delta2(x)

  std::array<El::Matrix<El::BigFloat>, 2> bilinear_bases;

  bilinear_bases[0]
    = sample_bilinear_basis(bilinear_basis[0], sample_points, sample_scalings);

  // The \sqrt(x) factors can be accounted for by replacing the
  // scale factors s_k with x_k s_k.
  std::vector<El::BigFloat> scaled_samples;
  for(size_t ii = 0; ii < sample_points.size(); ++ii)
    {
      scaled_samples.emplace_back(sample_points[ii] * sample_scalings[ii]);
    }
  bilinear_bases[1]
    = sample_bilinear_basis(bilinear_basis[1], sample_points, scaled_samples);
  return bilinear_bases;
}