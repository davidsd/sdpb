#include "../../../../sdp_convert.hxx"
#include "../../../../SDP.hxx"

void set_block_sizes(
  const std::vector<Dual_Constraint_Group> &dualConstraintGroups, SDP &sdp)
{
  for(auto &g : dualConstraintGroups)
    {
      sdp.dimensions.push_back(g.dim);
      sdp.degrees.push_back(g.degree);

      sdp.schur_block_sizes.push_back((g.dim * (g.dim + 1) / 2)
                                      * (g.degree + 1));

      // sdp.bilinear_bases is the concatenation of the g.bilinear_bases.
      // The matrix Y is a BlockDiagonalMatrix built from the
      // concatenation of the blocks for each individual
      // Dual_Constraint_Group.  sdp.blocks[j] = {b1, b2, ... } contains
      // the indices for the blocks of Y corresponding to the j-th
      // group.
      if(g.bilinear_bases.size() != 2)
        {
          throw std::runtime_error("Wrong number of elements in "
                                   "dualConstraintGroups::bilinear_bases.  "
                                   "Expected 2 but found: "
                                   + std::to_string(g.bilinear_bases.size()));
        }
      for(auto &b : g.bilinear_bases)
        {
          sdp.psd_matrix_block_sizes.push_back(b.Height() * g.dim);
          sdp.bilinear_pairing_block_sizes.push_back(b.Width() * g.dim);
        }
    }
}
