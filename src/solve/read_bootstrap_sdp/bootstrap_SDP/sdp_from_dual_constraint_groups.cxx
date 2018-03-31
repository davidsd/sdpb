#include "Dual_Constraint_Group.hxx"
#include "../../SDP.hxx"

// Collect a bunch of Dual_Constraint_Group's and a dual objective
// function into an SDP.

SDP sdp_from_dual_constraint_groups(
  const Vector &dualObjective, const Real &objectiveConst,
  const std::vector<Dual_Constraint_Group> &dualConstraintGroups)
{
  SDP sdp;
  sdp.dualObjective = dualObjective;
  sdp.objectiveConst = objectiveConst;

  for(std::vector<Dual_Constraint_Group>::const_iterator g
      = dualConstraintGroups.begin();
      g != dualConstraintGroups.end(); g++)
    {
      sdp.dimensions.push_back(g->dim);
      sdp.degrees.push_back(g->degree);

      // sdp.primalObjective is the concatenation of the
      // g.constraintConstants
      sdp.primalObjective.insert(sdp.primalObjective.end(),
                                 g->constraintConstants.begin(),
                                 g->constraintConstants.end());
    }
  sdp.FreeVarMatrix
    = Matrix(sdp.primalObjective.size(), sdp.dualObjective.size());

  int p = 0;
  // Each g corresponds to an index 0 <= j < J (not used explicitly here)
  for(std::vector<Dual_Constraint_Group>::const_iterator g
      = dualConstraintGroups.begin();
      g != dualConstraintGroups.end(); g++)
    {
      // sdp.bilinearBases is the concatenation of the g.bilinearBases.
      // The matrix Y is a BlockDiagonalMatrix built from the
      // concatenation of the blocks for each individual
      // Dual_Constraint_Group.  sdp.blocks[j] = {b1, b2, ... } contains
      // the indices for the blocks of Y corresponding to the j-th
      // group.
      std::vector<int> blocks;
      for(std::vector<Matrix>::const_iterator b = g->bilinearBases.begin();
          b != g->bilinearBases.end(); b++)
        {
          // Ensure that each bilinearBasis is sampled the correct number
          // of times
          assert(b->cols == g->degree + 1);
          blocks.push_back(sdp.bilinearBases.size());
          sdp.bilinearBases.push_back(*b);
        }
      sdp.blocks.push_back(blocks);

      // sdp.FreeVarMatrix is the block-wise concatenation of the
      // g.constraintMatrix's
      for(int k = 0; k < g->constraintMatrix.rows; k++, p++)
        for(int n = 0; n < g->constraintMatrix.cols; n++)
          sdp.FreeVarMatrix.elt(p, n) = g->constraintMatrix.elt(k, n);
    }
  assert(p == static_cast<int>(sdp.primalObjective.size()));

  sdp.initializeConstraintIndices();
  return sdp;
}
