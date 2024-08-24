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

#include "pmp/Polynomial_Vector_Matrix.hxx"
#include "sdpb_util/boost_serialization.hxx"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>

class Dual_Constraint_Group
{
public:
  size_t block_index{};
  size_t dim{};
  size_t num_points{};

  // constraint_matrix = B, a P'xN Matrix
  El::Matrix<El::BigFloat> constraint_matrix;

  // constraint_constants = c, a vector of length P'
  std::vector<El::BigFloat> constraint_constants;

  // bilinear_bases is a vector of Matrices encoding the v_{b,k}
  // entering the constraint matrices A_p, as described
  // above. `bilinear_bases' here has the structure of
  // `bilinear_bases[j]' above for some fixed j.
  std::array<El::Matrix<El::BigFloat>, 2> bilinear_bases;

  // Vector of length P', having the same structure as c.
  // Vector c and each column of matrix B are elementwise multiplied
  // by preconditioning_values:
  // c_{j,r,s,k} -> c_{j,r,s,k} * pv_{j,r,s,k}
  // B_{j,r,s,k},n -> B_{j,r,s,k},n * pv_{j,r,s,k}
  //
  // The goal of preconditioning is to make the final functional (c-B.y) flat,
  // thus reducing condition numbers etc.
  //
  // Currently, preconditioning_values{j,r,s,k} are calculated as (1+x_k)^{d_r + d_s},
  // where x_k are sample points, and d_i is a vector read from "preconditioningPowersVector" in pmp.json.
  //
  // preconditioning_values allows for more fine-grained adjustments
  // than prefactor/sample scalings (which have the same value for all PMP matrix elements)
  std::vector<El::BigFloat> preconditioning_values;

  Dual_Constraint_Group() = default;
  Dual_Constraint_Group(const size_t &Block_index,
                        const Polynomial_Vector_Matrix &m);
};
