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
