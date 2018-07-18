#include "../set_stream_precision.hxx"
#include "Dual_Constraint_Group.hxx"
#include "write_vector.hxx"

void write_primal_objective(
  const boost::filesystem::path &output_dir,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups)
{
  size_t block_index(0);

  for(auto &group : dual_constraint_groups)
    {
      assert(static_cast<size_t>(group.constraintMatrix.Height())
             == group.constraintConstants.size());

      boost::filesystem::ofstream output_stream(
        output_dir / ("primal_objective." + std::to_string(block_index)));
      set_stream_precision(output_stream);
      write_vector(output_stream, group.constraintConstants);
      ++block_index;
    }
}
