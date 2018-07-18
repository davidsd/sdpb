//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include "Block_Matrix.hxx"
#include "Block_Vector.hxx"
#include "../Polynomial_Vector_Matrix.hxx"
#include "Index_Tuple.hxx"
#include "ostream.hxx"

#include <boost/filesystem.hpp>

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

class SDP
{
public:
  // bilinearBases is a vector of Matrices encoding the v_{b,k}
  // that enter the constraint matrices. Specifically, v_{b,k} are
  // the columns of bilinearBases[b],
  //
  // bilinearBases[b].elt(m,k) = (v_{b,k})_m  (0 <= b < bMax,
  //                                           0 <= k <= d_j,
  //                                           0 <= m <= delta_b)
  //
  std::vector<El::Matrix<El::BigFloat>> bilinear_bases_local;
  std::vector<El::DistMatrix<El::BigFloat>> bilinear_bases_dist;

  // free_var_matrix = B, a PxN matrix
  Block_Matrix free_var_matrix;

  // c, a vector of length P used with primal_objective
  Block_Vector primal_objective_c;

  // b, a vector of length N used with dual_objective
  // It is duplicated amongst all the blocks
  El::DistMatrix<El::BigFloat> dual_objective_b;

  // objectiveConst = f
  El::BigFloat objective_const;

  // dimensions[j] = m_j  (0 <= j < J)
  std::vector<size_t> dimensions;

  // degrees[j] = d_j  (0 <= j < J)
  std::vector<size_t> degrees;

  SDP(const boost::filesystem::path &sdp_directory);

  std::vector<size_t> schur_block_sizes,
    // Dimensions of the blocks of X,Y
    // psd_matrix_block_sizes[b] = (delta_b+1)*m_j = length(v_{b,*})*m_j
    // (0 <= b < bMax)
    psd_matrix_block_sizes,
    // Dimensions of the bilinear pairing matrices U and V
    // bilinear_pairing_block_sizes[b] = (d_j + 1)*m_j
    // (0 <= b < bMax)
    bilinear_pairing_block_sizes;

  std::vector<size_t> block_indices;
  El::Grid grid;

  size_t total_psd_rows() const
  {
    return std::accumulate(psd_matrix_block_sizes.begin(),
                           psd_matrix_block_sizes.end(), 0);
  }

  std::vector<size_t> schur_offsets() const
  {
    size_t current_offset(0);
    std::vector<size_t> result;
    for(auto &size : schur_block_sizes)
      {
        result.push_back(current_offset);
        current_offset += size;
      }
    result.push_back(current_offset);
    return result;
  }
};
