#include "Dual_Constraint_Group.hxx"
#include "write_vector.hxx"
#include "../set_stream_precision.hxx"

void write_primal_objective_c(
  const boost::filesystem::path &output_dir,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups)
{
  size_t block_index(0);

  for(auto &group : dual_constraint_groups)
    {
      assert(static_cast<size_t>(group.constraint_matrix.Height())
             == group.constraint_constants.size());

      boost::filesystem::ofstream output_stream(
        output_dir / ("primal_objective_c." + std::to_string(block_index)));
      set_stream_precision(output_stream);
      write_vector(output_stream, group.constraint_constants);
      ++block_index;
    }
}
