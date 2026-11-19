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
#include "sdpb_util/assert.hxx"

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

// The n x K matrix with entries v_m(x_k) = sqrt(w_k) q_m(x_k), 0 <= m < n,
// where q_0, q_1, ... are orthonormal for sum_k w_k delta(x - x_k), made
// without polynomial coefficients by the Stieltjes (Lanczos) recurrence
//
//   beta_{m+1} v_{m+1} = x v_m - alpha_m v_m - beta_m v_{m-1},
//
// applied to the vectors themselves, with full reorthogonalization (twice)
// so that the rows stay orthonormal to the working precision.  Each row is
// sqrt(w) times a polynomial of degree m in x with a positive leading
// coefficient, as for the Cholesky-based bilinear_basis().  Used with
// pmp2sdp --encodingPoints (two_set_encoding.cxx): there the encoding points
// range over many decades and the monomial moments' Cholesky factorization
// loses most of the working precision.
namespace
{
  El::Matrix<El::BigFloat>
  sample_orthonormal_by_recurrence(const size_t n,
                                   const std::vector<El::BigFloat> &x,
                                   const std::vector<El::BigFloat> &w)
  {
    const size_t K = x.size();
    ASSERT_EQUAL(w.size(), K);
    ASSERT(n <= K, DEBUG_STRING(n), DEBUG_STRING(K));
    El::Matrix<El::BigFloat> v(n, K);
    if(n == 0)
      return v;

    const auto dot = [&](const size_t i, const std::vector<El::BigFloat> &r) {
      El::BigFloat s = 0;
      for(size_t k = 0; k < K; ++k)
        s += v.CRef(i, k) * r[k];
      return s;
    };
    const auto normalize_into = [&](const size_t row,
                                    const std::vector<El::BigFloat> &r) {
      El::BigFloat norm2 = 0;
      for(size_t k = 0; k < K; ++k)
        norm2 += r[k] * r[k];
      ASSERT(norm2 > 0, "bilinear basis by recurrence: zero vector at row ",
             row, DEBUG_STRING(n), DEBUG_STRING(K));
      const El::BigFloat inv_norm = 1 / El::Sqrt(norm2);
      for(size_t k = 0; k < K; ++k)
        v.Set(row, k, r[k] * inv_norm);
    };

    std::vector<El::BigFloat> r(K);
    for(size_t k = 0; k < K; ++k)
      {
        ASSERT(w[k] >= 0, DEBUG_STRING(k), DEBUG_STRING(w[k]));
        r[k] = El::Sqrt(w[k]);
      }
    normalize_into(0, r);

    for(size_t m = 0; m + 1 < n; ++m)
      {
        for(size_t k = 0; k < K; ++k)
          r[k] = x[k] * v.CRef(m, k);
        // Orthogonalize x v_m against v_0 .. v_m.  In exact arithmetic only
        // the coefficients of v_m (alpha_m) and v_{m-1} (beta_m) are
        // nonzero: this is the three-term recurrence.  The other terms keep
        // the rows orthonormal in finite precision.  Classical
        // Gram-Schmidt, twice.
        for(int pass = 0; pass < 2; ++pass)
          for(size_t j = 0; j <= m; ++j)
            {
              const El::BigFloat c = dot(j, r);
              for(size_t k = 0; k < K; ++k)
                r[k] -= c * v.CRef(j, k);
            }
        normalize_into(m + 1, r);
      }
    return v;
  }
}

std::array<El::Matrix<El::BigFloat>, 2> sample_bilinear_basis_by_recurrence(
  const std::array<size_t, 2> &basis_sizes,
  const std::vector<El::BigFloat> &sample_points,
  const std::vector<El::BigFloat> &sample_scalings)
{
  // Y_1: sqrt(s_k) q_m(x_k), q_m orthonormal for sum_k s_k delta(x - x_k);
  // Y_2: sqrt(x_k s_k) q_m(x_k), orthonormal for sum_k x_k s_k delta(x - x_k).
  std::vector<El::BigFloat> scaled_samples;
  for(size_t ii = 0; ii < sample_points.size(); ++ii)
    scaled_samples.emplace_back(sample_points[ii] * sample_scalings[ii]);
  return {sample_orthonormal_by_recurrence(basis_sizes[0], sample_points,
                                           sample_scalings),
          sample_orthonormal_by_recurrence(basis_sizes[1], sample_points,
                                           scaled_samples)};
}
