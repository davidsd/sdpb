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
      if(g.bilinearBases_elemental.size() != 2)
        {
          throw std::runtime_error(
            "Wrong number of elements in dualConstraintGroups::bilinearBases. "
            " "
            "Expected 2 but found:"
            + std::to_string(g.bilinearBases_elemental.size()));
        }
      sdp.blocks.push_back({sdp.bilinear_bases_elemental_local.size(),
                            sdp.bilinear_bases_elemental_local.size() + 1});
      for(auto &b : g.bilinearBases_elemental)
        {
          // Ensure that each bilinearBasis is sampled the correct number
          // of times
          assert(static_cast<size_t>(b.Width()) == g.degree + 1);
          sdp.bilinear_bases_elemental_local.push_back(b);
          sdp.bilinear_bases_elemental_dist.emplace_back(b.Height(),
                                                         b.Width());

          auto dist(sdp.bilinear_bases_elemental_dist.rbegin());
          for(int64_t row = 0; row < dist->LocalHeight(); ++row)
            {
              El::Int global_row(dist->GlobalRow(row));
              for(int64_t column = 0; column < dist->LocalWidth(); ++column)
                {
                  El::Int global_column(dist->GlobalCol(column));
                  dist->SetLocal(row, column,
                                 b.Get(global_row, global_column));
                }
            }
        }
    }

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
        if(block->GlobalCol(0) == 0)
          {
            for(size_t row = 0; row < local_height; ++row)
              {
                size_t global_row(block->GlobalRow(row));
                block->SetLocal(
                  row, 0, group->constraintConstants_elemental.at(global_row));
              }
          }
      }
      {
        sdp.free_var_matrix_elemental.blocks.emplace_back(
          block_size, sdp.dual_objective_b_elemental.Height());
        auto block(sdp.free_var_matrix_elemental.blocks.rbegin());
        int64_t local_height(block->LocalHeight()),
          local_width(block->LocalWidth());

        for(int64_t row = 0; row < local_height; ++row)
          {
            El::Int global_row(block->GlobalRow(row));
            if(global_row < group->constraintMatrix_elemental.Height())
              {
                for(int64_t column = 0; column < local_width; ++column)
                  {
                    El::Int global_column(block->GlobalCol(column));
                    if(global_column
                       < group->constraintMatrix_elemental.Width())
                      {
                        block->SetLocal(row, column,
                                        group->constraintMatrix_elemental.Get(
                                          global_row, global_column));
                      }
                  }
              }
          }
      }
      ++group;
    }
}
