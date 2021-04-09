#include "Dual_Constraint_Group.hxx"

#include <iostream>

void write_bilinear_bases(std::ostream &output_stream,
                          const Dual_Constraint_Group &group)
{
  output_stream << "  \"bilinear_bases_even\":\n  [\n";
  for(auto basis(group.bilinear_bases.begin());
      basis != group.bilinear_bases.end(); ++basis)
    {
      if(basis != group.bilinear_bases.begin())
        {
          output_stream << ",\n  \"bilinear_bases_odd\":\n  [\n";
        }
      // Ensure that each bilinearBasis is sampled the correct number
      // of times
      assert(static_cast<size_t>(basis->Width()) == group.degree + 1);
      for(int64_t row = 0; row < basis->Height(); ++row)
        {
          if(row != 0)
            {
              output_stream << ",\n";
            }
          output_stream << "    [\n";
          for(int64_t column = 0; column < basis->Width(); ++column)
            {
              if(column != 0)
                {
                  output_stream << ",\n";
                }
              output_stream << "      \"" << (*basis)(row, column) << "\"";
            }
          output_stream << "\n    ]";
        }
      output_stream << "\n  ]";
    }
  output_stream << ",\n";
}
