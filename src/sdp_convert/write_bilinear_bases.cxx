#include "Dual_Constraint_Group.hxx"
#include "../set_stream_precision.hxx"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <vector>

void write_bilinear_bases(
  const boost::filesystem::path &output_dir, const int &rank,
  const std::vector<size_t> &indices,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups)
{
  const boost::filesystem::path output_path(
    output_dir / ("bilinear_bases_" + std::to_string(rank) + ".json"));
  boost::filesystem::ofstream output_stream(output_path);
  set_stream_precision(output_stream);
  output_stream << "[\n";
  for(size_t block(0); block != dual_constraint_groups.size(); ++block)
    {
      if(block != 0)
        {
          output_stream << ",\n";
        }
      output_stream << "  {\n    \"block_index\": " << indices.at(block)
                    << ",\n    \"bilinear_bases_even\":\n    [\n";
      for(auto basis(dual_constraint_groups[block].bilinear_bases.begin());
          basis != dual_constraint_groups[block].bilinear_bases.end(); ++basis)
        {
          if(basis != dual_constraint_groups[block].bilinear_bases.begin())
            {
              output_stream << ",\n    \"bilinear_bases_odd\":\n    [\n";
            }
          // Ensure that each bilinearBasis is sampled the correct number
          // of times
          assert(static_cast<size_t>(basis->Width())
                 == dual_constraint_groups[block].degree + 1);
          for(int64_t row = 0; row < basis->Height(); ++row)
            {
              if(row != 0)
                {
                  output_stream << ",\n";
                }
              output_stream << "      [\n";
              for(int64_t column = 0; column < basis->Width(); ++column)
                {
                  if(column != 0)
                    {
                      output_stream << ",\n";
                    }
                  output_stream << "        \"" << (*basis)(row, column)
                                << "\"";
                }
              output_stream << "\n      ]";
            }
          output_stream << "\n    ]";
        }
      output_stream << "\n  }";
    }
  output_stream << "\n]\n";
  if(!output_stream.good())
    {
      throw std::runtime_error("Error when writing to: "
                               + output_path.string());
    }
}
