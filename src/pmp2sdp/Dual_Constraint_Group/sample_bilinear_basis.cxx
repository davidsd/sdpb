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

#include "pmp/Polynomial.hxx"

// Old algorithm
El::Matrix<El::BigFloat>
sample_bilinear_basis(const int maxDegree, const int numSamples,
                      const Polynomial_Vector &bilinearBasis,
                      const std::vector<El::BigFloat> &samplePoints,
                      const std::vector<El::BigFloat> &sampleScalings)
{
  El::Matrix<El::BigFloat> b(maxDegree + 1, numSamples);
  for(int k = 0; k < numSamples; k++)
    {
      El::BigFloat x(samplePoints.at(k)), scale(Sqrt(sampleScalings.at(k)));
      for(int i = 0; i <= maxDegree; i++)
        {
          b.Set(i, k, scale * bilinearBasis.at(i)(x));
        }
    }
  return b;
}

// New algorithm, if bilinearBasis not specified
std::array<El::Matrix<El::BigFloat>, 2>
sample_bilinear_basis(const std::vector<El::BigFloat> &sample_points,
                      const std::vector<El::BigFloat> &sample_scalings)
{
  std::array<El::Matrix<El::BigFloat>, 2> bilinear_bases;
  const size_t num_points = sample_points.size();
  ASSERT_EQUAL(sample_points.size(), sample_scalings.size());
  const size_t degree = num_points - 1;
  const std::array<size_t, 2> delta = {degree / 2, (degree + 1) / 2 - 1};

  for(const size_t index : {0, 1})
    {
      // Matrix of monomials
      // index = 0: set m[i,k] = scalings[k] * x[k]^i
      // index = 1: set m[i,k] = scalings[k] * x[k]^{i+1/2}
      El::Matrix<El::BigFloat> monomials(delta[index] + 1, num_points);
      for(int k = 0; k < monomials.Width(); ++k)
        {
          const auto &scaling = sample_scalings.at(k);
          const auto &x = sample_points.at(k);
          El::BigFloat x_pow = index == 0 ? 1 : El::Sqrt(x);
          for(int i = 0; i < monomials.Height(); ++i)
            {
              monomials.Set(i, k, El::Sqrt(scaling) * x_pow);
              x_pow *= x;
            }
        }

      // SVD decomposition of monomials
      El::Matrix<El::BigFloat> U, S, V;
      El::SVDCtrl<El::BigFloat> ctrl;
      ctrl.bidiagSVDCtrl.approach = El::COMPACT_SVD;
      // TODO set other ctrl parameters, as in compute_lamda.cxx?
      El::SVD(monomials, U, S, V);
      El::Transpose(V, bilinear_bases[index]);
    }

  return bilinear_bases;
}
