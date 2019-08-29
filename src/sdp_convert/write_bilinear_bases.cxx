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
    output_dir / ("bilinear_bases." + std::to_string(rank)));
  boost::filesystem::ofstream output_stream(output_path);
  set_stream_precision(output_stream);
  output_stream << dual_constraint_groups.size() << "\n";

  for(auto &group : dual_constraint_groups)
    {
      for(auto &basis : group.bilinear_bases)
        {
          // Ensure that each bilinearBasis is sampled the correct number
          // of times
          assert(static_cast<size_t>(basis.Width()) == group.degree + 1);
          output_stream << basis.Height() << " " << basis.Width() << "\n";
          for(int64_t row = 0; row < basis.Height(); ++row)
            for(int64_t column = 0; column < basis.Width(); ++column)
              {
                output_stream << basis(row, column) << "\n";
              }
        }
    }
  if(!output_stream.good())
    {
      throw std::runtime_error("Error when writing to: "
                               + output_path.string());
    }
}
