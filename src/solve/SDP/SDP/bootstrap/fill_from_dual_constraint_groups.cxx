#include "Dual_Constraint_Group.hxx"
#include "../../../SDP.hxx"

// Collect a bunch of Dual_Constraint_Group's and a dual objective
// function into an SDP.

void fill_from_dual_constraint_groups(
  const std::vector<Dual_Constraint_Group> &dualConstraintGroups, SDP &sdp)
{
  // Each g corresponds to an index 0 <= j < J (not used explicitly here)
  for(auto &g : dualConstraintGroups)
    {
      sdp.dimensions.push_back(g.dim);
      sdp.degrees.push_back(g.degree);

      // sdp.bilinear_bases is the concatenation of the g.bilinearBases.
      // The matrix Y is a BlockDiagonalMatrix built from the
      // concatenation of the blocks for each individual
      // Dual_Constraint_Group.  sdp.blocks[j] = {b1, b2, ... } contains
      // the indices for the blocks of Y corresponding to the j-th
      // group.
      std::vector<int> blocks;
      for(auto &b : g.bilinearBases_elemental)
        {
          // Ensure that each bilinearBasis is sampled the correct number
          // of times
          assert(static_cast<size_t>(b.Width()) == g.degree + 1);
          blocks.push_back(sdp.bilinear_bases_elemental.size());
          sdp.bilinear_bases_elemental.push_back(b);

          sdp.primal_objective_c_elemental.blocks.emplace_back(
            b.Height() * g.dim, 1);
          auto last_block(sdp.primal_objective_c_elemental.blocks.rbegin());
          size_t local_height(last_block->LocalHeight());
          El::Int row_min(last_block->GlobalRow(0));
          for(size_t hh = 0; hh < local_height; ++hh)
            {
              last_block->SetLocal(
                hh, 0, g.constraintConstants_elemental.at(row_min + hh));
            }
        }

      for(auto &b : g.bilinearBases)
        {
          // Ensure that each bilinearBasis is sampled the correct number
          // of times
          assert(b.cols == g.degree + 1);
          sdp.bilinear_bases.push_back(b);
        }
      sdp.blocks.push_back(blocks);

      // sdp.primal_objective is the concatenation of the
      // g.constraintConstants
      sdp.primal_objective_c.insert(sdp.primal_objective_c.end(),
                                    g.constraintConstants.begin(),
                                    g.constraintConstants.end());
    }

  sdp.free_var_matrix = Matrix(sdp.primal_objective_c.size(),
                               sdp.dual_objective_b_elemental.Height());
  size_t row = 0;
  for(auto &g : dualConstraintGroups)
    {
      // sdp.free_var_matrix is the block-wise concatenation of the
      // g.constraintMatrix's
      for(size_t k = 0; k < g.constraintMatrix.rows; ++k, ++row)
        for(size_t n = 0; n < g.constraintMatrix.cols; ++n)
          {
            sdp.free_var_matrix.elt(row, n) = g.constraintMatrix.elt(k, n);
          }
    }
  assert(row == sdp.primal_objective_c.size());

  sdp.initialize_constraint_indices();
}
