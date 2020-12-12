#include "Dual_Constraint_Group.hxx"
#include "../set_stream_precision.hxx"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <vector>

void write_bilinear_bases(
  const boost::filesystem::path &output_dir, const int &rank,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups)
{
  const boost::filesystem::path output_path(
    output_dir / ("bilinear_bases_" + std::to_string(rank) + ".json"));
  boost::filesystem::ofstream output_stream(output_path);
  set_stream_precision(output_stream);
  output_stream << "{\n  \"dual_constraint_groups\": [\n";
  for(auto group(dual_constraint_groups.begin());
      group != dual_constraint_groups.end(); ++group)
    {
      if(group != dual_constraint_groups.begin())
        {
          output_stream << ",\n";
        }
      output_stream << "    {\n      \"bilinear_bases\": [\n";
      for(auto basis(group->bilinear_bases.begin());
          basis != group->bilinear_bases.end(); ++basis)
        {
          if(basis != group->bilinear_bases.begin())
            {
              output_stream << ",\n";
            }
          // Ensure that each bilinearBasis is sampled the correct number
          // of times
          assert(static_cast<size_t>(basis->Width()) == group->degree + 1);
          output_stream << "        [\n";
          for(int64_t row = 0; row < basis->Height(); ++row)
            {
              if(row != 0)
                {
                  output_stream << ",\n";
                }
              output_stream << "          [\n";
              for(int64_t column = 0; column < basis->Width(); ++column)
                {
                  if(column != 0)
                    {
                      output_stream << ",\n";
                    }
                  output_stream << "            \"" << (*basis)(row, column)
                                << "\"";
                }
              output_stream << "\n          ]";
            }
          output_stream << "\n        ]";
        }
      output_stream << "\n      ]\n    }";
    }
  output_stream << "\n  ]\n}\n";
  if(!output_stream.good())
    {
      throw std::runtime_error("Error when writing to: "
                               + output_path.string());
    }
}
