#include "../Dual_Constraint_Group.hxx"
#include "../../../../SDP.hxx"

void set_block_sizes_and_bilinear(
  const std::vector<Dual_Constraint_Group> &dualConstraintGroups, SDP &sdp);
void assign_blocks(
  const std::vector<Dual_Constraint_Group> &dualConstraintGroups, SDP &sdp);

// Collect a bunch of Dual_Constraint_Group's and a dual objective
// function into an SDP.

void fill_from_dual_constraint_groups(
  const std::vector<Dual_Constraint_Group> &dualConstraintGroups, SDP &sdp)
{
  set_block_sizes_and_bilinear(dualConstraintGroups, sdp);
  assign_blocks(dualConstraintGroups, sdp);
}
