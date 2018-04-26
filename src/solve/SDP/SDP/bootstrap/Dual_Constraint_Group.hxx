#pragma once

// Dual_Constraint_Group represents a set of constraints of the form
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
// Dual_Constraint_Groups labeled by 0<=j<J.
//
// Dual_Constraint_Group's are currently only used as an intermediate
// data structure between the polynomial matrices defining an PMP and a
// full SDP.  By directly combining Dual_Constraint_Groups into an SDP
// using sdpFromDual_Constraint_Groups, it is possible to define slightly
// more general optimization problems than those produced by
// bootstrapSDP.  Perhaps this level of generality will be useful in
// the future.
//

#include "../../../Matrix.hxx"

#include <El.hpp>

class Dual_Constraint_Group
{
public:
  size_t dim;
  size_t degree;

  // constraintMatrix = B, a P'xN Matrix
  Matrix constraintMatrix;
  El::Matrix<El::BigFloat> constraintMatrix_elemental;

  // constraintConstants = c, a vector of length P'
  Vector constraintConstants;
  std::vector<El::BigFloat> constraintConstants_elemental;

  // bilinearBases is a vector of Matrices encoding the v_{b,k}
  // entering the constraint matrices A_p, as described
  // above. `bilinearBases' here has the structure of
  // `bilinearBases[j]' above for some fixed j.
  std::vector<Matrix> bilinearBases;
  std::vector<El::Matrix<El::BigFloat>> bilinearBases_elemental;
};
