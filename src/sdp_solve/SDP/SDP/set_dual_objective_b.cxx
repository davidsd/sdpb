#include <El.hpp>

void set_dual_objective_b(const std::vector<El::BigFloat> &objective_input,
                          const El::Grid &grid,
                          El::DistMatrix<El::BigFloat> &dual_objective_b)
{
  dual_objective_b.SetGrid(grid);
  dual_objective_b.Resize(objective_input.size(), 1);
  if(dual_objective_b.GlobalCol(0) == 0)
    {
      size_t local_height(dual_objective_b.LocalHeight());
      for(size_t row = 0; row < local_height; ++row)
        {
          size_t global_row(dual_objective_b.GlobalRow(row));
          dual_objective_b.SetLocal(row, 0, objective_input.at(global_row));
        }
    }
}
