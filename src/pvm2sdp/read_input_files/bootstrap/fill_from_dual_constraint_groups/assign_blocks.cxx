#include "../../../../sdp_convert.hxx"
#include "../../../../SDP.hxx"

void assign_blocks(
  const std::vector<Dual_Constraint_Group> &dualConstraintGroups, SDP &sdp)
{
  sdp.primal_objective_c.blocks.reserve(sdp.block_indices.size());
  sdp.free_var_matrix.blocks.reserve(sdp.block_indices.size());
  
  for(auto &block_index : sdp.block_indices)
    {
      const auto &block_size(sdp.schur_block_sizes.at(block_index));
      const auto &group(dualConstraintGroups.at(block_index));

      assert(group.constraintConstants.size() == block_size);
      assert(static_cast<size_t>(group.constraintMatrix.Height())
             == block_size);
      {
        sdp.primal_objective_c.blocks.emplace_back(block_size, 1, sdp.grid);
        auto &block(sdp.primal_objective_c.blocks.back());
        size_t local_height(block.LocalHeight());
        if(block.GlobalCol(0) == 0)
          {
            for(size_t row = 0; row < local_height; ++row)
              {
                size_t global_row(block.GlobalRow(row));
                block.SetLocal(row, 0,
                               group.constraintConstants.at(global_row));
              }
          }
      }
      {
        sdp.free_var_matrix.blocks.emplace_back(
          block_size, sdp.dual_objective_b.Height(), sdp.grid);
        auto &block(sdp.free_var_matrix.blocks.back());
        int64_t local_height(block.LocalHeight()),
          local_width(block.LocalWidth());

        for(int64_t row = 0; row < local_height; ++row)
          {
            El::Int global_row(block.GlobalRow(row));
            if(global_row < group.constraintMatrix.Height())
              {
                for(int64_t column = 0; column < local_width; ++column)
                  {
                    El::Int global_column(block.GlobalCol(column));
                    if(global_column < group.constraintMatrix.Width())
                      {
                        block.SetLocal(row, column,
                                       group.constraintMatrix.Get(
                                         global_row, global_column));
                      }
                  }
              }
          }
      }
    }
}
