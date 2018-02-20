//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#include <vector>
#include "SDP.h"

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
Matrix sampleBilinearBasis(const int maxDegree,
                           const int numSamples,
                           const vector<Polynomial> &bilinearBasis,
                           const vector<Real> &samplePoints,
                           const vector<Real> &sampleScalings) {
  Matrix b(maxDegree + 1, numSamples);
  for (int k = 0; k < numSamples; k++) {
    Real x     = samplePoints[k];
    Real scale = sqrt(sampleScalings[k]);
    for (int i = 0; i <= maxDegree; i++)
      b.elt(i, k) = scale*bilinearBasis[i](x);
  }
  return b;
}

// Convert a PolynomialVectorMatrix to a DualConstraint group by
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
DualConstraintGroup
dualConstraintGroupFromPolVecMat(const PolynomialVectorMatrix &m) {
  DualConstraintGroup g;

  assert(m.rows == m.cols);
  g.dim    = m.rows;
  g.degree = m.degree();

  int numSamples     = g.degree + 1;
  int numConstraints = numSamples * g.dim * (g.dim + 1)/2;
  int vectorDim      = m.elt(0, 0).size();

  // Form the constraintMatrix B and constraintConstants c from the
  // polynomials (1,y) . \vec P^{rs}(x)

  // The first element of each vector \vec P^{rs}(x) multiplies the constant 1
  g.constraintConstants = Vector(numConstraints);
  // The rest multiply decision variables y
  g.constraintMatrix    = Matrix(numConstraints, vectorDim - 1);

  // Populate B and c by sampling the polynomial matrix
  int p = 0;
  for (int c = 0; c < g.dim; c++) {
    for (int r = 0; r <= c; r++) {
      for (int k = 0; k < numSamples; k++) {
        Real x     = m.samplePoints[k];
        Real scale = m.sampleScalings[k];

        g.constraintConstants[p] = scale*m.elt(r, c)[0](x);
        for (int n = 1; n < vectorDim; n++)
          g.constraintMatrix.elt(p, n-1) = -scale*m.elt(r, c)[n](x);

        p++;
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
  int delta1 = g.degree/2;
  g.bilinearBases.push_back(sampleBilinearBasis(delta1, numSamples,
                                                m.bilinearBasis,
                                                m.samplePoints,
                                                m.sampleScalings));
  int delta2 = (g.degree - 1)/2;
  // a degree-0 PolynomialVectorMatrix only needs one block
  if (delta2 >= 0)
    // The \sqrt(x) factors can be accounted for by replacing the
    // scale factors s_k with x_k s_k.
    g.bilinearBases
      .push_back(sampleBilinearBasis(delta2, numSamples,
                                     m.bilinearBasis,
                                     m.samplePoints,
                                     multiplyVectors(m.samplePoints,
                                                     m.sampleScalings)));

  return g;
}

// Collect a bunch of DualConstraintGroup's and a dual objective
// function into an SDP.
SDP sdpFromDualConstraintGroups(const Vector &dualObjective,
                                const Real &objectiveConst,
                                const vector<DualConstraintGroup> &dualConstraintGroups) {
  SDP sdp;
  sdp.dualObjective  = dualObjective;
  sdp.objectiveConst = objectiveConst;

  for (vector<DualConstraintGroup>::const_iterator g = dualConstraintGroups.begin();
       g != dualConstraintGroups.end();
       g++) {
    sdp.dimensions.push_back(g->dim);
    sdp.degrees.push_back(g->degree);

    // sdp.primalObjective is the concatenation of the
    // g.constraintConstants
    sdp.primalObjective.insert(sdp.primalObjective.end(),
                               g->constraintConstants.begin(),
                               g->constraintConstants.end());
  }
  sdp.FreeVarMatrix = Matrix(sdp.primalObjective.size(), sdp.dualObjective.size());

  int p = 0;
  // Each g corresponds to an index 0 <= j < J (not used explicitly here)
  for (vector<DualConstraintGroup>::const_iterator g = dualConstraintGroups.begin();
       g != dualConstraintGroups.end();
       g++) {
    // sdp.bilinearBases is the concatenation of the g.bilinearBases.
    // The matrix Y is a BlockDiagonalMatrix built from the
    // concatenation of the blocks for each individual
    // DualConstraintGroup.  sdp.blocks[j] = {b1, b2, ... } contains
    // the indices for the blocks of Y corresponding to the j-th
    // group.
    vector<int> blocks;
    for (vector<Matrix>::const_iterator b = g->bilinearBases.begin();
         b != g->bilinearBases.end();
         b++) {
      // Ensure that each bilinearBasis is sampled the correct number
      // of times
      assert(b->cols == g->degree + 1);
      blocks.push_back(sdp.bilinearBases.size());
      sdp.bilinearBases.push_back(*b);
    }
    sdp.blocks.push_back(blocks);

    // sdp.FreeVarMatrix is the block-wise concatenation of the
    // g.constraintMatrix's
    for (int k = 0; k < g->constraintMatrix.rows; k++, p++)
      for (int n = 0; n < g->constraintMatrix.cols; n++)
        sdp.FreeVarMatrix.elt(p, n) = g->constraintMatrix.elt(k, n);
  }
  assert(p == static_cast<int>(sdp.primalObjective.size()));

  sdp.initializeConstraintIndices();
  return sdp;
}

// Form an SDP from an affineObjective and a list of
// PolynomialVectorMatrices.  affineObjective is a vector of the form
// (f, b), where the objective function will be
//
//   f + y.b
//
// with y the dual decision variables.  polVectorMatrices are
// described in SDP.h.
//
SDP bootstrapSDP(const Vector &affineObjective,
                 const vector<PolynomialVectorMatrix> &polVectorMatrices) {
  // Convert polVectorMatrices into DualConstraintGroup's
  vector<DualConstraintGroup> dualConstraintGroups;
  for (vector<PolynomialVectorMatrix>::const_iterator m = polVectorMatrices.begin();
       m != polVectorMatrices.end(); m++)
    dualConstraintGroups.push_back(dualConstraintGroupFromPolVecMat(*m));

  // Split affineObjective into objectiveConst f and dualObjective b
  Real objectiveConst = affineObjective[0];
  Vector dualObjective = affineObjective;
  dualObjective.erase(dualObjective.begin());

  return sdpFromDualConstraintGroups(dualObjective, objectiveConst, dualConstraintGroups);
}

