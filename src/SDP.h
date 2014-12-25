//=======================================================================
// Copyright 2014 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#ifndef SDPB_SDP_H_
#define SDPB_SDP_H_

#include <algorithm>
#include <vector>
#include <iostream>
#include <ostream>
#include "types.h"
#include "util.h"
#include "Vector.h"
#include "Matrix.h"
#include "Polynomial.h"

using std::vector;
using std::ostream;

// The constraint matrices A_* are labeled by tuples
// (j_1,r_1,s_1,k_1), (j_2,r_2,s_2,k_2), ...
// The corresponding IndexTuples are 
// IndexTuple(p=1,r=r_1,s=s_1,k=k_1), IndexTuple(p=2,r=r_2,s=s_2,k=k_2), ...
class IndexTuple {
 public:
  int p;
  int r;
  int s;
  int k;
  IndexTuple(int p, int r, int s, int k): p(p), r(r), s(s), k(k) {}
  IndexTuple() {}
};

class SDP {
 public:
  vector<Matrix> bilinearBases;
  Matrix FreeVarMatrix;
  Vector primalObjective;
  Vector dualObjective;
  Real objectiveConst;
  vector<int> dimensions;
  vector<int> degrees;
  vector<vector<int> > blocks;
  vector<vector<IndexTuple> > constraintIndices;

  vector<int> psdMatrixBlockDims() const {
    vector<int> dims;
    for (unsigned int j = 0; j < dimensions.size(); j++)
      for (vector<int>::const_iterator b = blocks[j].begin();
           b != blocks[j].end(); b++)
        dims.push_back(bilinearBases[*b].rows * dimensions[j]);
    return dims;
  }

  vector<int> bilinearPairingBlockDims() const {
    vector<int> dims;
    for (unsigned int j = 0; j < dimensions.size(); j++)
      for (vector<int>::const_iterator b = blocks[j].begin();
           b != blocks[j].end(); b++)
        dims.push_back(bilinearBases[*b].cols * dimensions[j]);
    return dims;
  }

  vector<int> schurBlockDims() const {
    vector<int> dims;
    for (unsigned int j = 0; j < dimensions.size(); j++)
      dims.push_back(constraintIndices[j].size());
    return dims;
  }

  void initializeConstraintIndices() {
    int p = 0;
    for (unsigned int j = 0; j < dimensions.size(); j++) {
      constraintIndices.push_back(vector<IndexTuple>(0));

      for (int s = 0; s < dimensions[j]; s++) {
        for (int r = 0; r <= s; r++) {
          for (int k = 0; k <= degrees[j]; k++) {
            constraintIndices[j].push_back(IndexTuple(p, r, s, k));
            p++;
          }
        }
      }
    }
    assert(p == static_cast<int>(primalObjective.size()));
  }

  friend ostream& operator<<(ostream& os, const SDP& sdp) {
    os << "SDP(bilinearBases = " << sdp.bilinearBases
       << ", FreeVarMatrix = " << sdp.FreeVarMatrix
       << ", primalObjective = " << sdp.primalObjective
       << ", dualObjective = " << sdp.dualObjective
       << ", dimensions = " << sdp.dimensions
       << ", degrees = " << sdp.degrees
       << ", blocks = " << sdp.blocks
       << ")";
    return os;
  }
};

class SampledMatrixPolynomial {
 public:
  int dim;
  int degree;
  Matrix constraintMatrix;
  Vector constraintConstants;
  vector<Matrix> bilinearBases;
};

class PolynomialVectorMatrix {
 public:
  int rows;
  int cols;
  vector<vector<Polynomial> > elements;
  vector<Real> samplePoints;
  vector<Real> sampleScalings;
  vector<Polynomial> bilinearBasis;

  inline const vector<Polynomial>& elt(const int r, const int c) const {
    return elements[r + c*rows];
  }

  inline vector<Polynomial>& elt(const int r, const int c) {
    return elements[r + c*rows];
  }

  int degree() const {
    int d = 0;
    for (vector<vector<Polynomial> >::const_iterator e = elements.begin();
         e != elements.end(); e++)
      for (vector<Polynomial>::const_iterator p = e->begin();
           p != e->end(); p++)
        d = max(p->degree(), d);
    return d;
  }
};

SampledMatrixPolynomial samplePolynomialVectorMatrix(const PolynomialVectorMatrix &m);

SDP bootstrapSDP(const Vector &objective,
                 const Real &objectiveConst,
                 const vector<SampledMatrixPolynomial> &sampledMatrixPols);

SDP bootstrapPolynomialSDP(const Vector &affineObjective,
                           const vector<PolynomialVectorMatrix> &polVectorMatrices);

#endif  // SDPB_SDP_H_
