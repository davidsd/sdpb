#include "Dual_Constraint_Group.hxx"
#include "write_vector.hxx"

void write_blocks(
  const boost::filesystem::path &output_dir, const int &rank,
  const std::vector<size_t> &indices,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups)
{
  const boost::filesystem::path output_path(
    output_dir / ("blocks_" + std::to_string(rank) + ".json"));
  boost::filesystem::ofstream output_stream(output_path);
  output_stream << "{\n  \"blocks\": [\n";
  for(size_t block(0); block != indices.size(); ++block)
    {
      if(block != 0)
        {
          output_stream << ",\n";
        }
      output_stream << "    {\n      \"index\": " << indices.at(block)
                    << ",\n      \"dim\": "
                    << dual_constraint_groups.at(block).dim
                    << ",\n      \"num_points\": "
                    << (dual_constraint_groups.at(block).degree + 1)
                    << "\n    }";
    }
  output_stream << "\n  ]\n}\n";
  if(!output_stream.good())
    {
      throw std::runtime_error("Error when writing to: "
                               + output_path.string());
    }
}
