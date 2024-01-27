//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#include "../Dual_Constraint_Group.hxx"
#include "sdpb_util/assert.hxx"

void write_block_json(std::ostream &output_stream,
                      const Dual_Constraint_Group &group);

El::Matrix<El::BigFloat>
sample_bilinear_basis(const int maxDegree, const int numSamples,
                      const Polynomial_Vector &bilinearBasis,
                      const std::vector<El::BigFloat> &samplePoints,
                      const std::vector<El::BigFloat> &sampleScalings);

// Construct a Dual_Constraint_Group from a Polynomial_Vector_Matrix by
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
Dual_Constraint_Group::Dual_Constraint_Group(
  const size_t &Block_index, const Polynomial_Vector_Matrix &m)
    : block_index(Block_index),
      dim(m.polynomials.Height()),
      num_points(m.sample_points.size())
{
  const auto &polys = m.polynomials;
  ASSERT(polys.Height() == polys.Width());

  const size_t degree(num_points - 1),
    numConstraints(num_points * dim * (dim + 1) / 2),
    vectorDim(polys(0, 0).size());

  // Form the constraint_matrix B and constraint_constants c from the
  // polynomials (1,y) . \vec P^{rs}(x)

  // The first element of each vector \vec P^{rs}(x) multiplies the constant 1
  constraint_constants.resize(numConstraints);
  // The rest multiply decision variables y
  constraint_matrix.Resize(numConstraints, vectorDim - 1);

  // Populate B and c by sampling the polynomial matrix
  int p = 0;
  for(size_t c = 0; c < dim; c++)
    {
      for(size_t r = 0; r <= c; r++)
        {
          for(size_t k = 0; k < num_points; k++)
            {
              El::BigFloat x(m.sample_points.at(k));
              El::BigFloat scale(m.sample_scalings.at(k));
              constraint_constants[p] = scale * polys(r, c)[0](x);
              for(size_t n = 1; n < vectorDim; ++n)
                {
                  constraint_matrix.Set(p, n - 1, -scale * polys(r, c)[n](x));
                }
              ++p;
            }
        }
    }

  // The matrix Y has two blocks Y_1, Y_2.  The bilinear_bases for the
  // constraint matrices A_p are given by sampling the following
  // vectors for each block:
  //
  //   Y_1: {q_0(x), ..., q_delta1(x)}
  //   Y_2: {\sqrt(x) q_0(x), ..., \sqrt(x) q_delta2(x)
  //
  const size_t delta1(degree / 2);
  bilinear_bases[0] = sample_bilinear_basis(
    delta1, num_points, m.bilinear_basis, m.sample_points, m.sample_scalings);

  // For degree==0, the second block will have zero size.
  const size_t delta2((degree + 1) / 2 - 1);
  // The \sqrt(x) factors can be accounted for by replacing the
  // scale factors s_k with x_k s_k.
  std::vector<El::BigFloat> scaled_samples;
  for(size_t ii = 0; ii < m.sample_points.size(); ++ii)
    {
      scaled_samples.emplace_back(m.sample_points[ii] * m.sample_scalings[ii]);
    }
  bilinear_bases[1] = sample_bilinear_basis(
    delta2, num_points, m.bilinear_basis, m.sample_points, scaled_samples);
}
