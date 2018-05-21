//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#include "../Dual_Constraint_Group.hxx"
#include "../../../../SDP.hxx"

El::Matrix<El::BigFloat>
sample_bilinear_basis(const int maxDegree, const int numSamples,
                      const std::vector<Polynomial> &bilinearBasis,
                      const std::vector<El::BigFloat> &samplePoints,
                      const std::vector<El::BigFloat> &sampleScalings);

// Convert a Polynomial_Vector_Matrix to a DualConstraint group by
// sampling the matrix at the appropriate number of points, as
// described in SDP.h:
//
//   (1,y) . M(x) is positive semidefinite
//
// is equivalent to
//
//   Tr(A_p Y) + (B y)_p = c_p
//
// for tuples p = (r,s,k).
//
Dual_Constraint_Group
dual_constraint_group_from_pol_vec_mat(const Polynomial_Vector_Matrix &m)
{
  Dual_Constraint_Group g;

  assert(m.rows == m.cols);
  g.dim = m.rows;
  g.degree = m.degree();

  size_t numSamples = g.degree + 1;
  size_t numConstraints = numSamples * g.dim * (g.dim + 1) / 2;
  size_t vectorDim = m.elt(0, 0).size();

  // Form the constraintMatrix B and constraintConstants c from the
  // polynomials (1,y) . \vec P^{rs}(x)

  // The first element of each vector \vec P^{rs}(x) multiplies the constant 1
  g.constraintConstants = Vector(numConstraints);
  g.constraintConstants_elemental.resize(numConstraints);
  // The rest multiply decision variables y
  g.constraintMatrix = Matrix(numConstraints, vectorDim - 1);
  g.constraintMatrix_elemental.Resize(numConstraints, vectorDim - 1);

  // Populate B and c by sampling the polynomial matrix
  int p = 0;
  for(size_t c = 0; c < g.dim; c++)
    {
      for(size_t r = 0; r <= c; r++)
        {
          for(size_t k = 0; k < numSamples; k++)
            {
              Real x = m.sample_points[k];
              Real scale = m.sample_scalings[k];

              g.constraintConstants[p] = scale * m.elt(r, c)[0](x);
              for(size_t n = 1; n < vectorDim; ++n)
                {
                  g.constraintMatrix.elt(p, n - 1)
                    = -scale * m.elt(r, c)[n](x);
                }

              {
                El::BigFloat x = m.sample_points_elemental[k];
                El::BigFloat scale = m.sample_scalings_elemental[k];
                g.constraintConstants_elemental[p] = scale * m.elt(r, c)[0](x);
                for(size_t n = 1; n < vectorDim; ++n)
                  {
                    g.constraintMatrix_elemental.Set(
                      p, n - 1, -scale * m.elt(r, c)[n](x));
                  }
              }
              ++p;
            }
        }
    }

  // The matrix Y has two blocks Y_1, Y_2.  The bilinearBases for the
  // constraint matrices A_p are given by sampling the following
  // vectors for each block:
  //
  //   Y_1: {q_0(x), ..., q_delta1(x)}
  //   Y_2: {\sqrt(x) q_0(x), ..., \sqrt(x) q_delta2(x)
  //
  int delta1 = g.degree / 2;
  g.bilinearBases_elemental.push_back(sample_bilinear_basis(
    delta1, numSamples, m.bilinear_basis, m.sample_points_elemental,
    m.sample_scalings_elemental));

  int delta2 = (g.degree - 1) / 2;
  // a degree-0 Polynomial_Vector_Matrix only needs one block
  if(delta2 >= 0)
    {
      // The \sqrt(x) factors can be accounted for by replacing the
      // scale factors s_k with x_k s_k.
      std::vector<El::BigFloat> scaled_samples;
      for(size_t ii = 0; ii < m.sample_points_elemental.size(); ++ii)
        {
          scaled_samples.emplace_back(m.sample_points_elemental[ii]
                                      * m.sample_scalings_elemental[ii]);
        }
      g.bilinearBases_elemental.push_back(
        sample_bilinear_basis(delta2, numSamples, m.bilinear_basis,
                              m.sample_points_elemental, scaled_samples));
    }
  return g;
}
