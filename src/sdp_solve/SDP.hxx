//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include "Block_Info.hxx"
#include "Block_Matrix.hxx"
#include "Block_Vector.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <filesystem>
#include <optional>

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

struct SDP
{
  // bilinear_bases is a vector of Matrices encoding the v_{b,k}
  // that enter the constraint matrices. Specifically, v_{b,k} are
  // the columns of bilinear_bases[b],
  //
  // bilinear_bases[b].elt(m,k) = (v_{b,k})_m  (0 <= b < bMax,
  //                                           0 <= k <= d_j,
  //                                           0 <= m <= delta_b)

  std::vector<El::DistMatrix<El::BigFloat>> bilinear_bases;

  // Blocks of bilinear bases precomputed for X_inv and Y
  std::vector<El::DistMatrix<El::BigFloat>> bases_blocks;

  // free_var_matrix = B, a PxN matrix
  Block_Matrix free_var_matrix;

  // c, a vector of length P used with primal_objective
  Block_Vector primal_objective_c;

  // b, a vector of length N used with dual_objective
  // It is duplicated amongst all the blocks
  El::DistMatrix<El::BigFloat> dual_objective_b;

  // objectiveConst = f
  El::BigFloat objective_const;

  // Vector of length N+1, initialized only on rank=0
  // Normalization condition n.z = 1 is used to translate
  // SDP_Solver.y (eq. 2.2 in SDPB Manual) to vector z (eq. 3.1)
  // Not initialized if there was no normalization in PMP
  std::optional<std::vector<El::BigFloat>> normalization;

  SDP(const std::filesystem::path &sdp_path, const Block_Info &block_info,
      const El::Grid &grid, Timers &timers);
  SDP(const El::BigFloat &objective_const,
      const std::vector<std::vector<El::BigFloat>> &primal_objective_c_input,
      const std::vector<El::Matrix<El::BigFloat>> &free_var_input,
      const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &yp_to_y_star,
      const El::DistMatrix<El::BigFloat, El::STAR, El::STAR>
        &dual_objective_b_star,
      const std::vector<El::BigFloat> &normalization,
      const El::BigFloat &primal_c_scale, const Block_Info &block_info,
      const El::Grid &grid);

private:
  void validate(const Block_Info &block_info) const noexcept(false);
};
