#include "Dual_Constraint_Group.hxx"
#include "write_vector.hxx"

void write_primal_objective_c(std::ostream &output_stream,
                              const Dual_Constraint_Group &group)
{
  assert(static_cast<size_t>(group.constraint_matrix.Height())
         == group.constraint_constants.size());

  output_stream << "  \"c\":\n";
  write_vector(output_stream, group.constraint_constants, "  ");
  output_stream << ",\n";
}
