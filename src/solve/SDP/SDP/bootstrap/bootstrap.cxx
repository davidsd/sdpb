//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#include "Dual_Constraint_Group.hxx"
#include "../../../SDP.hxx"

void fill_from_dual_constraint_groups(
  const std::vector<Dual_Constraint_Group> &dualConstraintGroups, SDP &sdp);

Dual_Constraint_Group
dual_constraint_group_from_pol_vec_mat(const Polynomial_Vector_Matrix &m);

// Form an SDP from an affineObjective and a list of
// PolynomialVectorMatrices.  affineObjective is a vector of the form
// (f, b), where the objective function will be
//
//   f + y.b
//
// with y the dual decision variables.  polVectorMatrices are
// described in SDP.h.
//
void bootstrap(const std::vector<El::BigFloat> &affine_objective_elemental,
               const std::vector<Polynomial_Vector_Matrix> &polVectorMatrices,
               SDP &sdp)
{
  // Convert polVectorMatrices into Dual_Constraint_Group's
  std::vector<Dual_Constraint_Group> dualConstraintGroups;
  for(auto &m : polVectorMatrices)
    {
      dualConstraintGroups.push_back(
        dual_constraint_group_from_pol_vec_mat(m));
    }

  // Split affine_objective into objectiveConst f and dualObjective b
  auto affine(affine_objective_elemental.begin());
  assert(affine != affine_objective_elemental.end());
  sdp.objective_const_elemental = *affine;
  ++affine;

  sdp.dual_objective_b_elemental.Resize(affine_objective_elemental.size() - 1,
                                        1);
  size_t local_height(sdp.dual_objective_b_elemental.LocalHeight());

  if(sdp.dual_objective_b_elemental.GlobalCol(0) == 0)
    {
      for(size_t row = 0; row < local_height; ++row)
        {
          size_t global_row(sdp.dual_objective_b_elemental.GlobalRow(row));
          sdp.dual_objective_b_elemental.SetLocal(
            row, 0, affine_objective_elemental[global_row + 1]);
        }
    }
  fill_from_dual_constraint_groups(dualConstraintGroups, sdp);
}
