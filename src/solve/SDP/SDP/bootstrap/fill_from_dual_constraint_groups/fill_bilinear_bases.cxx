#include "../Dual_Constraint_Group.hxx"
#include "../../../../SDP.hxx"

void fill_bilinear_bases(
  const std::vector<Dual_Constraint_Group> &dualConstraintGroups, SDP &sdp)
{
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
      if(g.bilinearBases.size() != 2)
        {
          throw std::runtime_error("Wrong number of elements in "
                                   "dualConstraintGroups::bilinearBases.  "
                                   "Expected 2 but found:"
                                   + std::to_string(g.bilinearBases.size()));
        }
      sdp.blocks.push_back({sdp.bilinear_bases_local.size(),
                            sdp.bilinear_bases_local.size() + 1});
      for(auto &b : g.bilinearBases)
        {
          // Ensure that each bilinearBasis is sampled the correct number
          // of times
          assert(static_cast<size_t>(b.Width()) == g.degree + 1);
          sdp.bilinear_bases_local.push_back(b);
          sdp.bilinear_bases_dist.emplace_back(b.Height(), b.Width());

          auto dist(sdp.bilinear_bases_dist.rbegin());
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
}
