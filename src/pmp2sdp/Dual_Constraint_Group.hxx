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

  // Two-set encoding (pmp2sdp --encodingPoints): decimal digits the rows
  // lose in the step from the approximation points to the encoding points
  // (estimate, see two_set_encoding.cxx); not written to the SDP.
  double encoding_digits_lost = 0;

  Dual_Constraint_Group() = default;
  Dual_Constraint_Group(const size_t &Block_index,
                        const Polynomial_Vector_Matrix &m);
};

// The second step of the two-set encoding (two_set_encoding.cxx): rewrite a
// block written at its approximation points (rows sigma_k R(x_k)) at the
// given encoding points and scalings, with bilinear bases by recurrence.
// Returns the estimated decimal digits lost.
double reencode_dual_constraint_group(
  Dual_Constraint_Group &group,
  const std::vector<El::BigFloat> &approximation_points,
  const std::vector<El::BigFloat> &approximation_scalings,
  const std::vector<El::BigFloat> &encoding_points,
  const std::vector<El::BigFloat> &encoding_scalings,
  const std::array<size_t, 2> &basis_sizes);

// The scale of the blocks of degree 0 under the two-set encoding
// (scale_degree_0_blocks.cxx): the rows of B and c of every degree-0 block
// (one sample point) are multiplied by the power of 2 that brings their
// largest |entry|, over b^2 (b its bilinear basis), nearest to the median,
// over the re-encoded blocks of all ranks, of their largest |entry| (an
// equivalent constraint; the basis is unchanged).  Its sample scalings are
// multiplied by the same power, so pmp_info.json stays consistent with the
// rows.  The counts and exponents are over all ranks.
struct Degree_0_Scaling
{
  El::Int num_degree_0 = 0;
  El::Int num_reencoded = 0;
  // Blocks multiplied by 2^e with e != 0, and the range of e (if
  // num_degree_0 > 0 and num_reencoded > 0)
  El::Int num_scaled = 0;
  El::Int min_exponent = 0, max_exponent = 0;
  // Median of log2(largest |entry|) over the re-encoded blocks
  double median_log2 = 0;
};
// groups[i] must be the block of matrices[i]; comm: the ranks that hold the
// blocks (collective).
Degree_0_Scaling
scale_degree_0_blocks(std::vector<Dual_Constraint_Group> &groups,
                      std::vector<Polynomial_Vector_Matrix> &matrices,
                      const El::mpi::Comm &comm = El::mpi::COMM_WORLD);
// The median of the given values (the even case: the mean of the two middle
// ones).
double median_log2_max_abs(std::vector<double> values);
