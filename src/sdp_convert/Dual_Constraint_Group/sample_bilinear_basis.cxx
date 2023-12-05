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

#include "sdpb_util/Polynomial.hxx"

El::Matrix<El::BigFloat>
sample_bilinear_basis(const int maxDegree, const int numSamples,
                      const std::vector<Polynomial> &bilinearBasis,
                      const std::vector<El::BigFloat> &samplePoints,
                      const std::vector<El::BigFloat> &sampleScalings)
{
  El::Matrix<El::BigFloat> b(maxDegree + 1, numSamples);
  for(int k = 0; k < numSamples; k++)
    {
      El::BigFloat x(samplePoints[k]), scale(Sqrt(sampleScalings[k]));
      for(int i = 0; i <= maxDegree; i++)
        {
          b.Set(i, k, scale * bilinearBasis[i](x));
        }
    }
  return b;
}
