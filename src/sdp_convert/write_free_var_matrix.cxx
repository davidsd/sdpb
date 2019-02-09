#include "Dual_Constraint_Group.hxx"
#include "write_vector.hxx"
#include "../set_stream_precision.hxx"

void write_free_var_matrix(
  const boost::filesystem::path &output_dir,
  const size_t &dual_objectives_b_size,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups)
{
  size_t block_index(0);

  for(auto &group : dual_constraint_groups)
    {
      size_t block_size(group.constraint_matrix.Height());

      boost::filesystem::ofstream output_stream(
        output_dir / ("free_var_matrix." + std::to_string(block_index)));
      set_stream_precision(output_stream);

      output_stream << block_size << " " << dual_objectives_b_size << "\n";
      for(size_t row = 0; row < block_size; ++row)
        for(size_t column = 0; column < dual_objectives_b_size; ++column)
          {
            output_stream << group.constraint_matrix(row, column) << "\n";
          }

      ++block_index;
    }
}
