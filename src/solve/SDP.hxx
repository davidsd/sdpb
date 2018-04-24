//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include "Matrix.hxx"
#include "Polynomial_Vector_Matrix.hxx"
#include "Index_Tuple.hxx"

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
  std::vector<Matrix> bilinear_bases;

  // FreeVarMatrix = B, a PxN matrix
  Matrix free_var_matrix;

  // c, a vector of length P used with primal_objective
  Vector primal_objective_c;

  // b, a vector of length N used with dual_objective
  Vector dual_objective_b;

  // objectiveConst = f
  Real objective_const;

  // dimensions[j] = m_j  (0 <= j < J)
  std::vector<int> dimensions;

  El::BigFloat objective_const_elemental;
  // degrees[j] = d_j  (0 <= j < J)
  std::vector<int> degrees;

  // blocks gives the 1-to-many mapping
  //
  // blocks[j] = {b_{j1}, b_{j2}, ...}  (0 <= j < J)
  //
  // entering the constraint matrices A_p.  blocks[j1] and blocks[j2]
  // are disjoint unless j1==j2.
  //
  std::vector<std::vector<int>> blocks;

  // constraintIndices gives the 1-to-many mapping described above
  //
  // constraintIndices[j] = { Index_Tuple(p,r,s,k)
  //                          for 0 <= s < m_j,
  //                              0 <= r <= s,
  //                              0 <= k <= d_j,
  //                          with p the overall constraint index }
  //
  // This allows us to loop through the constraints A_p associated to
  // each j.
  //
  std::vector<std::vector<Index_Tuple>> constraint_indices;

  SDP(const std::vector<boost::filesystem::path> &sdp_files);

  // create the mapping j -> { Index_Tuple(p,r,s,k) } described above
  //
  void initialize_constraint_indices()
  {
    int p = 0;
    for(unsigned int j = 0; j < dimensions.size(); j++)
      {
        constraint_indices.emplace_back(0);

        for(int s = 0; s < dimensions[j]; s++)
          for(int r = 0; r <= s; r++)
            for(int k = 0; k <= degrees[j]; k++)
              {
                constraint_indices[j].emplace_back(p, r, s, k);
                p++;
              }
      }
    assert(p == static_cast<int>(primal_objective_c.size()));
  }

  // Dimensions of the blocks of X,Y (0 <= b < bMax)
  //
  // psdMatrixBlockDims()[b] = (delta_b+1)*m_j = length(v_{b,*})*m_j
  //
  std::vector<int> psd_matrix_block_dims() const
  {
    std::vector<int> dims;
    for(size_t j = 0; j < dimensions.size(); j++)
      for(auto &b : blocks[j])
        {
          dims.push_back(bilinear_bases[b].rows * dimensions[j]);
        }
    return dims;
  }

  // Dimensions of the bilinear pairing matrices U^(b) and V^(b) (0 <= b < bMax)
  //
  // bilinearPairingBlockDims()[b] = (d_j + 1)*m_j
  //
  std::vector<int> bilinear_pairing_block_dims() const
  {
    std::vector<int> dims;
    for(size_t j = 0; j < dimensions.size(); j++)
      for(auto &b : blocks[j])
        {
          dims.push_back(bilinear_bases[b].cols * dimensions[j]);
        }
    return dims;
  }

  // Dimensions of the blocks S^(j) of the Schur complement matrix:
  //
  // schurBlockDims()[j] = (d_j+1)*m_j*(m_j+1)/2
  //                     = length(constraintIndices[j])
  //
  std::vector<int> schur_block_dims() const
  {
    std::vector<int> dims;
    for(unsigned int j = 0; j < dimensions.size(); j++)
      dims.push_back(constraint_indices[j].size());
    return dims;
  }

  // Print an SDP, for debugging purposes
  friend std::ostream &operator<<(std::ostream &os, const SDP &sdp)
  {
    os << "SDP(bilinearBases = " << sdp.bilinear_bases
       << ", FreeVarMatrix = " << sdp.free_var_matrix
       << ", primalObjective = " << sdp.primal_objective_c
       << ", dualObjective = " << sdp.dual_objective_b
       << ", dimensions = " << sdp.dimensions << ", degrees = " << sdp.degrees
       << ", blocks = " << sdp.blocks << ")";
    return os;
  }
};
