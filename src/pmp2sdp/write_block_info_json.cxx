#include "Dual_Constraint_Group.hxx"

#include <iostream>

void write_block_info_json(std::ostream &output_stream,
                           const Dual_Constraint_Group &group)
{
  output_stream << "{\n";
  output_stream << "  \"dim\": " << group.dim
                << ",\n  \"num_points\": " << (group.num_points);
  output_stream << "\n}";
}
