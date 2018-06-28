#include "../Dual_Constraint_Group.hxx"
#include "../../../../SDP.hxx"

void set_block_sizes(
  const std::vector<Dual_Constraint_Group> &dualConstraintGroups, SDP &sdp);
void set_bilinear(
  const std::vector<Dual_Constraint_Group> &dualConstraintGroups, SDP &sdp);
void compute_block_grid_mapping(const size_t &num_blocks,
                                std::list<El::Grid> &mapping);
void assign_blocks(
  const std::vector<Dual_Constraint_Group> &dualConstraintGroups,
  const std::list<El::Grid> &block_grid_mapping, SDP &sdp);

// Collect a bunch of Dual_Constraint_Group's and a dual objective
// function into an SDP.

void fill_from_dual_constraint_groups(
  const std::vector<El::BigFloat> &affine_objective,
  const std::vector<Dual_Constraint_Group> &dualConstraintGroups, SDP &sdp)
{
  set_block_sizes(dualConstraintGroups, sdp);
  compute_block_grid_mapping(sdp.dimensions.size(), sdp.block_grid_mapping);
  set_bilinear(dualConstraintGroups, sdp);

  // Split affine_objective into objectiveConst f and dualObjective b
  sdp.objective_const = affine_objective.at(0);

  // SetGrid() must happen before Resize()
  sdp.dual_objective_b.SetGrid(sdp.block_grid_mapping.front());
  sdp.dual_objective_b.Resize(affine_objective.size() - 1, 1);
  if(sdp.dual_objective_b.GlobalCol(0) == 0)
    {
      size_t local_height(sdp.dual_objective_b.LocalHeight());
      for(size_t row = 0; row < local_height; ++row)
        {
          size_t global_row(sdp.dual_objective_b.GlobalRow(row));
          sdp.dual_objective_b.SetLocal(row, 0,
                                        affine_objective[global_row + 1]);
        }
    }
  assign_blocks(dualConstraintGroups, sdp.block_grid_mapping, sdp);
}
