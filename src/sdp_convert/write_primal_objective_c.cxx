#include "Dual_Constraint_Group.hxx"
#include "write_vector.hxx"
#include "../set_stream_precision.hxx"

void write_primal_objective_c(
  const boost::filesystem::path &output_dir,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups)
{
  for(auto &group : dual_constraint_groups)
    {
      assert(static_cast<size_t>(group.constraint_matrix.Height())
             == group.constraint_constants.size());

      const boost::filesystem::path output_path(
        output_dir / ("primal_objective_c_" + std::to_string(group.block_index) + ".json"));
      boost::filesystem::ofstream output_stream(output_path);
      set_stream_precision(output_stream);
      write_vector(output_stream, group.constraint_constants, "");
      output_stream << "\n";
      if(!output_stream.good())
        {
          throw std::runtime_error("Error when writing to: "
                                   + output_path.string());
        }
    }
}
