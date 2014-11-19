//=======================================================================
// Copyright 2014 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#include <vector>
#include "SDP.h"

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

SampledMatrixPolynomial
samplePolynomialVectorMatrix(const PolynomialVectorMatrix &m) {
  SampledMatrixPolynomial s;

  assert(m.rows == m.cols);
  s.dim    = m.rows;
  s.degree = m.degree();

  int numSamples     = s.degree + 1;
  int numConstraints = numSamples * s.dim * (s.dim + 1)/2;
  int vectorDim      = m.elt(0, 0).size();


  // The first element of each vector multiplies the constant 1
  s.constraintConstants = Vector(numConstraints);
  // The rest multiply decision variables
  s.constraintMatrix    = Matrix(numConstraints, vectorDim - 1);

  int p = 0;
  for (int c = 0; c < s.dim; c++) {
    for (int r = c; r < s.dim; r++) {
      for (int k = 0; k < numSamples; k++) {
        Real x     = m.samplePoints[k];
        Real scale = m.sampleScalings[k];

        s.constraintConstants[p] = scale*m.elt(r, c)[0](x);
        for (int n = 1; n < vectorDim; n++)
          s.constraintMatrix.elt(p, n-1) = -scale*m.elt(r, c)[n](x);

        p++;
      }
    }
  }

  int delta1 = s.degree/2;
  s.bilinearBases.push_back(sampleBilinearBasis(delta1, numSamples,
                                                m.bilinearBasis,
                                                m.samplePoints,
                                                m.sampleScalings));
  int delta2 = (s.degree - 1)/2;
  if (delta2 >= 0)
    s.bilinearBases
      .push_back(sampleBilinearBasis(delta2, numSamples,
                                     m.bilinearBasis,
                                     m.samplePoints,
                                     multiplyVectors(m.samplePoints,
                                                     m.sampleScalings)));

  return s;
}

SDP bootstrapSDP(const Vector &objective,
                 const Real &objectiveConst,
                 const vector<SampledMatrixPolynomial> &sampledMatrixPols) {
  SDP sdp;
  sdp.dualObjective  = objective;
  sdp.objectiveConst = objectiveConst;

  for (vector<SampledMatrixPolynomial>::const_iterator s = sampledMatrixPols.begin();
       s != sampledMatrixPols.end();
       s++) {
    sdp.dimensions.push_back(s->dim);
    sdp.degrees.push_back(s->degree);
    sdp.primalObjective.insert(sdp.primalObjective.end(),
                               s->constraintConstants.begin(),
                               s->constraintConstants.end());
  }
  sdp.FreeVarMatrix = Matrix(sdp.primalObjective.size(), sdp.dualObjective.size());

  int p = 0;
  for (vector<SampledMatrixPolynomial>::const_iterator s = sampledMatrixPols.begin();
       s != sampledMatrixPols.end();
       s++) {
    vector<int> blocks;
    for (vector<Matrix>::const_iterator b = s->bilinearBases.begin();
         b != s->bilinearBases.end();
         b++) {
      assert(b->cols == s->degree + 1);
      blocks.push_back(sdp.bilinearBases.size());
      sdp.bilinearBases.push_back(*b);
    }
    sdp.blocks.push_back(blocks);

    for (int k = 0; k < s->constraintMatrix.rows; k++, p++)
      for (int n = 0; n < s->constraintMatrix.cols; n++)
        sdp.FreeVarMatrix.elt(p, n) = s->constraintMatrix.elt(k, n);
  }
  assert(p == static_cast<int>(sdp.primalObjective.size()));

  sdp.initializeConstraintIndices();
  return sdp;
}

SDP bootstrapPolynomialSDP(const Vector &affineObjective,
                           const vector<PolynomialVectorMatrix> &polVectorMatrices) {
  vector<SampledMatrixPolynomial> sampledMatrixPols;
  for (vector<PolynomialVectorMatrix>::const_iterator m = polVectorMatrices.begin();
       m != polVectorMatrices.end(); m++)
    sampledMatrixPols.push_back(samplePolynomialVectorMatrix(*m));

  Vector objective = affineObjective;
  objective.erase(objective.begin());
  Real objectiveConst = affineObjective[0];

  return bootstrapSDP(objective, objectiveConst, sampledMatrixPols);
}

