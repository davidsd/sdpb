#include "Dual_Constraint_Group.hxx"

#include <iostream>

void write_blocks(std::ostream &output_stream,
                  const Dual_Constraint_Group &group)
{
  output_stream << "  \"dim\": " << group.dim
                << ",\n  \"num_points\": " << (group.degree + 1) << ",\n";
}
