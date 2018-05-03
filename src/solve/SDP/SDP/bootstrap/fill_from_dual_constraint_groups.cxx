#include "Dual_Constraint_Group.hxx"
#include "../../../SDP.hxx"

// Collect a bunch of Dual_Constraint_Group's and a dual objective
// function into an SDP.

void fill_from_dual_constraint_groups(
  const std::vector<Dual_Constraint_Group> &dualConstraintGroups, SDP &sdp)
{
  // First compute blocks
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
      std::vector<size_t> blocks;
      for(auto &b : g.bilinearBases_elemental)
        {
          // Ensure that each bilinearBasis is sampled the correct number
          // of times
          assert(static_cast<size_t>(b.Width()) == g.degree + 1);
          blocks.push_back(sdp.bilinear_bases_elemental_local.size());
          sdp.bilinear_bases_elemental_local.push_back(b);
          sdp.bilinear_bases_elemental_dist.emplace_back(b.Height(),
                                                         b.Width());

          auto dist(sdp.bilinear_bases_elemental_dist.rbegin());
          El::Int row_offset(dist->GlobalRow(0)),
            column_offset(dist->GlobalCol(0));
          for(int64_t row = 0; row < dist->LocalHeight(); ++row)
            for(int64_t column = 0; column < dist->LocalWidth(); ++column)
              {
                dist->SetLocal(
                  row, column,
                  b.Get(row + row_offset, column + column_offset));
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

  // Then assign blocks
  auto group(dualConstraintGroups.begin());
  for(auto &block_size : sdp.schur_block_dims())
    {
      assert(group != dualConstraintGroups.end());
      assert(group->constraintConstants_elemental.size() == block_size);
      assert(static_cast<size_t>(group->constraintMatrix_elemental.Height())
             == block_size);
      {
        sdp.primal_objective_c_elemental.blocks.emplace_back(block_size, 1);
        auto block(sdp.primal_objective_c_elemental.blocks.rbegin());
        size_t local_height(block->LocalHeight());
        El::Int row_begin(block->GlobalRow(0));
        for(size_t hh = 0; hh < local_height; ++hh)
          {
            block->SetLocal(
              hh, 0, group->constraintConstants_elemental.at(row_begin + hh));
          }
      }
      {
        sdp.free_var_matrix_elemental.blocks.emplace_back(
          block_size, sdp.dual_objective_b_elemental.Height());
        auto block(sdp.free_var_matrix_elemental.blocks.rbegin());
        int64_t local_height(block->LocalHeight()),
          local_width(block->LocalWidth());
        El::Int row_begin(block->GlobalRow(0)),
          column_begin(block->GlobalCol(0));

        for(int64_t row = 0;
            row < local_height
            && row + row_begin < group->constraintMatrix_elemental.Height();
            ++row)
          for(int64_t column = 0;
              column < local_width
              && column + column_begin
                   < group->constraintMatrix_elemental.Width();
              ++column)
            {
              block->SetLocal(row, column,
                              group->constraintMatrix_elemental.Get(
                                row + row_begin, column + column_begin));
            }
      }
      ++group;
    }
}
