#include "Dual_Constraint_Group.hxx"

#include <iostream>

void write_free_var_matrix(std::ostream &output_stream,
                           const size_t &dual_objectives_b_size,
                           const Dual_Constraint_Group &group)
{
  size_t block_size(group.constraint_matrix.Height());

  output_stream << "  \"B\":\n  [\n";
  for(size_t row = 0; row < block_size; ++row)
    {
      if(row != 0)
        {
          output_stream << ",\n";
        }
      output_stream << "    [\n";
      for(size_t column = 0; column < dual_objectives_b_size; ++column)
        {
          if(column != 0)
            {
              output_stream << ",\n";
            }
          output_stream << "      \"" << group.constraint_matrix(row, column)
                        << "\"";
        }
      output_stream << "\n    ]";
    }
  output_stream << "\n  ]\n";
}
