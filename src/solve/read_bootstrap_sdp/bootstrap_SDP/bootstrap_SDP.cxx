//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#include "Dual_Constraint_Group.hxx"
#include "../../SDP.hxx"

Dual_Constraint_Group
dual_constraint_group_from_pol_vec_mat(const Polynomial_Vector_Matrix &m);

// Collect a bunch of Dual_Constraint_Group's and a dual objective
// function into an SDP.
SDP sdpFromDualConstraintGroups(
  const Vector &dualObjective, const Real &objectiveConst,
  const vector<Dual_Constraint_Group> &dualConstraintGroups)
{
  SDP sdp;
  sdp.dualObjective = dualObjective;
  sdp.objectiveConst = objectiveConst;

  for(vector<Dual_Constraint_Group>::const_iterator g
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
  for(vector<Dual_Constraint_Group>::const_iterator g
      = dualConstraintGroups.begin();
      g != dualConstraintGroups.end(); g++)
    {
      // sdp.bilinearBases is the concatenation of the g.bilinearBases.
      // The matrix Y is a BlockDiagonalMatrix built from the
      // concatenation of the blocks for each individual
      // Dual_Constraint_Group.  sdp.blocks[j] = {b1, b2, ... } contains
      // the indices for the blocks of Y corresponding to the j-th
      // group.
      vector<int> blocks;
      for(vector<Matrix>::const_iterator b = g->bilinearBases.begin();
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

// Form an SDP from an affineObjective and a list of
// PolynomialVectorMatrices.  affineObjective is a vector of the form
// (f, b), where the objective function will be
//
//   f + y.b
//
// with y the dual decision variables.  polVectorMatrices are
// described in SDP.h.
//
SDP bootstrap_SDP(const Vector &affineObjective,
                  const vector<Polynomial_Vector_Matrix> &polVectorMatrices)
{
  // Convert polVectorMatrices into Dual_Constraint_Group's
  vector<Dual_Constraint_Group> dualConstraintGroups;
  for(vector<Polynomial_Vector_Matrix>::const_iterator m
      = polVectorMatrices.begin();
      m != polVectorMatrices.end(); m++)
    dualConstraintGroups.push_back(dual_constraint_group_from_pol_vec_mat(*m));

  // Split affineObjective into objectiveConst f and dualObjective b
  Real objectiveConst = affineObjective[0];
  Vector dualObjective = affineObjective;
  dualObjective.erase(dualObjective.begin());

  return sdpFromDualConstraintGroups(dualObjective, objectiveConst,
                                     dualConstraintGroups);
}
