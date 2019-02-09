#include "../../../../sdp_convert.hxx"
#include "../../../../SDP.hxx"

void set_bilinear(
  const std::vector<Dual_Constraint_Group> &dualConstraintGroups, SDP &sdp)
{
  for(auto &g : dualConstraintGroups)
    {
      for(auto &b : g.bilinearBases)
        {
          // Ensure that each bilinearBasis is sampled the correct number
          // of times
          assert(static_cast<size_t>(b.Width()) == g.degree + 1);
          sdp.bilinear_bases_local.push_back(b);
          sdp.bilinear_bases_dist.emplace_back(b.Height(), b.Width(),
                                               sdp.grid);

          auto &dist(sdp.bilinear_bases_dist.back());
          for(int64_t row = 0; row < dist.LocalHeight(); ++row)
            {
              El::Int global_row(dist.GlobalRow(row));
              for(int64_t column = 0; column < dist.LocalWidth(); ++column)
                {
                  El::Int global_column(dist.GlobalCol(column));
                  dist.SetLocal(row, column, b.Get(global_row, global_column));
                }
            }
        }
    }
}
