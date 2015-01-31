//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
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
#include "Vector.h"
#include "Matrix.h"
#include "Polynomial.h"

using std::vector;
using std::ostream;

// The class SDP encodes a semidefinite program of the following form
//
// Dual: maximize f + b.y over y,Y such that
//                Tr(A_p Y) + (B y)_p = c_p  (0 <= p < P)
//                Y >= 0
// Primal: minimize f + c.x over x,X such that
//                X = \sum_p A_p x_p - C
//                B^T x = b
//                X >= 0
//
// where the data of the SDP has the following structure
//
// - b,y are vectors of length N
//
// - c,x are vectors of length P
//
// - f is a constant that we add to the objective functions for
//   convenience. It has no effect on the running of our algorithm.
//
// - B is a P x N matrix, called the free variable matrix
//
// - X and Y are block diagonal:
//
//   X = BlockDiagonal(X^(0), ..., X^(bMax-1))
//   Y = BlockDiagonal(Y^(0), ..., Y^(bMax-1))
//
//   Let us define Block_b(M) as BlockDiagonal(0,...,0,M,0,...,0),
//   where M is in the b-th block.  Then X and Y can be written
//
//   X = \sum_{0<=b<bMax} Block_b(X^(b))
//   Y = \sum_{0<=b<bMax} Block_b(Y^(b))
//
// - The constraints labeled by 0 <= p < P are in 1-to-1
//   correspondence with tuples
//   
//   p <-> (j,r,s,k) where 0 <= j < J,
//                         0 <= s < m_j,
//                         0 <= r <= s,
//                         0 <= k <= d_j,
//
//   We often interchange the index p with the tuple (j,r,s,k).  The
//   constraint matrices A_p are given by
//
//   A_(j,r,s,k) = \sum_{b \in blocks[j]}
//                     Block_b(v_{b,k} v_{b,k}^T \otimes E^{rs}),
//
//   where
//   - E^{rs} is the symmetrization of the m_j x m_j matrix with a
//     1 in entry (r,s) and zeros elsewhere.
//   - v_{b,k} is a vector of length (delta_b+1)
//   - \otimes denotes a tensor product
//   - each block above thus has dimension (delta_b+1)*m_j.
//   - blocks[j_1] and blocks[j_2] are disjoint if j_1 != j_2.  Thus,
//     the block index b determines a unique j, but not vice-versa.
//

// Instead of working with the 1 to 1 correspondence p <-> (j,r,s,k),
// it is convenient to collect tuples with the same index j. This
// gives a 1-to-many correspondence
//
// j <-> { list of IndexTuple(p,r,s,k) }
//
// IndexTuple is simply a named version of the 4-tuple (p,r,s,k).  A
// given IndexTuple uniquely specifies a constraint matrix A_p.
//
class IndexTuple {
 public:
  int p; // overall index of the constraint
  int r; // first index for E^{rs}
  int s; // second index for E^{rs}
  int k; // index for v_{b,k}
  IndexTuple(int p, int r, int s, int k): p(p), r(r), s(s), k(k) {}
  IndexTuple() {}
};

class SDP {
 public:
  // bilinearBases is a vector of Matrices encoding the v_{b,k}
  // that enter the constraint matrices. Specifically, v_{b,k} are
  // the columns of bilinearBases[b],
  //
  // bilinearBases[b].elt(m,k) = (v_{b,k})_m  (0 <= b < bMax,
  //                                           0 <= k <= d_j,
  //                                           0 <= m <= delta_b)
  //
  vector<Matrix> bilinearBases;

  // FreeVarMatrix = B, a PxN matrix
  Matrix FreeVarMatrix;

  // primalObjective = c, a vector of length P
  Vector primalObjective;

  // dualObjective = b, a vector of length N
  Vector dualObjective;

  // objectiveConst = f
  Real objectiveConst;

  // dimensions[j] = m_j  (0 <= j < J)
  vector<int> dimensions;

  // degrees[j] = d_j  (0 <= j < J)
  vector<int> degrees;

  // blocks gives the 1-to-many mapping
  //
  // blocks[j] = {b_{j1}, b_{j2}, ...}  (0 <= j < J)
  //
  // entering the constraint matrices A_p.  blocks[j1] and blocks[j2]
  // are disjoint unless j1==j2.
  //
  vector<vector<int> > blocks;

  // constraintIndices gives the 1-to-many mapping described above 
  //
  // constraintIndices[j] = { IndexTuple(p,r,s,k)
  //                          for 0 <= s < m_j,
  //                              0 <= r <= s,
  //                              0 <= k <= d_j,
  //                          with p the overall constraint index }
  //
  // This allows us to loop through the constraints A_p associated to
  // each j.
  //
  vector<vector<IndexTuple> > constraintIndices;

  // create the mapping j -> { IndexTuple(p,r,s,k) } described above
  //
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

  // Dimensions of the blocks of X,Y (0 <= b < bMax)
  //
  // psdMatrixBlockDims()[b] = (delta_b+1)*m_j = length(v_{b,*})*m_j
  //
  vector<int> psdMatrixBlockDims() const {
    vector<int> dims;
    for (unsigned int j = 0; j < dimensions.size(); j++)
      for (vector<int>::const_iterator b = blocks[j].begin();
           b != blocks[j].end(); b++)
        dims.push_back(bilinearBases[*b].rows * dimensions[j]);
    return dims;
  }

  // Dimensions of the bilinear pairing matrices U^(b) and V^(b) (0 <= b < bMax)
  //
  // bilinearPairingBlockDims()[b] = (d_j + 1)*m_j
  //
  vector<int> bilinearPairingBlockDims() const {
    vector<int> dims;
    for (unsigned int j = 0; j < dimensions.size(); j++)
      for (vector<int>::const_iterator b = blocks[j].begin();
           b != blocks[j].end(); b++)
        dims.push_back(bilinearBases[*b].cols * dimensions[j]);
    return dims;
  }

  // Dimensions of the blocks S^(j) of the Schur complement matrix:
  //
  // schurBlockDims()[j] = (d_j+1)*m_j*(m_j+1)/2
  //                     = length(constraintIndices[j])
  //
  vector<int> schurBlockDims() const {
    vector<int> dims;
    for (unsigned int j = 0; j < dimensions.size(); j++)
      dims.push_back(constraintIndices[j].size());
    return dims;
  }

  // Print an SDP, for debugging purposes
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

// DualConstraintGroup represents a set of constraints of the form
//
//   Tr(A_p Y) + (B y)_p = c_p
// 
// for a fixed j in the definition of SDP above. Here p corresponds to
//
//   p <-> (r,s,k) where 0 <= s < dim,
//                       0 <= r <= s,
//                       0 <= k <= degree
//
//   0 <= p < (degree+1)*dim*(dim+1)/2 = P'
//
// The constraints of a full SDP can be thought of as a collection of
// DualConstraintGroups labeled by 0<=j<J.
//
// DualConstraintGroup's are currently only used as an intermediate
// data structure between the polynomial matrices defining an PMP and a
// full SDP.  By directly combining DualConstraintGroups into an SDP
// using sdpFromDualConstraintGroups, it is possible to define slightly
// more general optimization problems than those produced by
// bootstrapSDP.  Perhaps this level of generality will be useful in
// the future.
//
class DualConstraintGroup {
 public:
  int dim;
  int degree;

  // constraintMatrix = B, a P'xN Matrix
  Matrix constraintMatrix;

  // constraintConstants = c, a vector of length P'
  Vector constraintConstants;

  // bilinearBases is a vector of Matrices encoding the v_{b,k}
  // entering the constraint matrices A_p, as described
  // above. `bilinearBases' here has the structure of
  // `bilinearBases[j]' above for some fixed j.
  vector<Matrix> bilinearBases;
};

// Let M(x) be a matrix whose entries are vectors of polynomials:
//
//   M(x) = ( \vec P^{00}(x) ... \vec P^{m0}(x) )
//          ( ...                               )
//          ( \vec P^{0n}(x) ... \vec P^{mn}(x) )
//
// where each vector has length N+1:
//
//   \vec P^{rs}(x) = (P^{rs}_{-1}(x), P^{rs}_0, ... , P^{rs}_{N-1}(x))
//
// Consider a vector y = (y_0, ..., y_{N-1}) of length N, and let
// (1,y) denote the vector of length N+1 whose components are 1,
// followed by the components of y.  As explained in the manual, the
// constraint
//
//   (1,y) . M(x) is positive semidefinite
//
// is equivalent to a DualConstraintGroup
//
//   Tr(A_p Y) + (B y)_p = c_p
//
// A PolynomialVectorMatrix contains the data needed to construct this
// DualConstraintGroup:  
//
class PolynomialVectorMatrix {
 public:
  int rows; // rows of M
  int cols; // cols of M

  // elements of M, in row-major order
  vector<vector<Polynomial> > elements;

  // A list of real numbers x_k (0 <= k <= degree(M)) at which to
  // sample M(x) to construct the v_{b,k}.
  vector<Real> samplePoints;

  // A list of real numbers s_k (0 <= k <= degree(M)) to scale M(x_k)
  // and the corresponding v_{b,k}.
  vector<Real> sampleScalings;

  // bilinearBasis[m] = q_m(x) (0 <= m <= degree/2), where q_m is a
  // polynomial with degree deg(q_m) = m.
  vector<Polynomial> bilinearBasis;

  inline const vector<Polynomial>& elt(const int r, const int c) const {
    return elements[r + c*rows];
  }

  inline vector<Polynomial>& elt(const int r, const int c) {
    return elements[r + c*rows];
  }

  // The maximal degree of any of the components P^{rs}_n(x).
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

SDP bootstrapSDP(const Vector &affineObjective,
                 const vector<PolynomialVectorMatrix> &polVectorMatrices);

#endif  // SDPB_SDP_H_
