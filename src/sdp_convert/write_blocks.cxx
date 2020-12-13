#include "Dual_Constraint_Group.hxx"
#include "write_vector.hxx"

void write_blocks(
  const boost::filesystem::path &output_dir,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups)
{
  for(auto &block : dual_constraint_groups)
    {
      const boost::filesystem::path output_path(
        output_dir
        / ("blocks_" + std::to_string(block.block_index) + ".json"));
      boost::filesystem::ofstream output_stream(output_path);
      output_stream << "{\n  \"dim\": " << block.dim
                    << ",\n  \"num_points\": " << (block.degree + 1)
                    << "\n}\n";
      if(!output_stream.good())
        {
          throw std::runtime_error("Error when writing to: "
                                   + output_path.string());
        }
    }
}
