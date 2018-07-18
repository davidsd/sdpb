#include "../set_stream_precision.hxx"

#include <boost/filesystem.hpp>
#include <vector>

void write_objective(const boost::filesystem::path &output_dir,
                     const std::vector<El::BigFloat> &objectives)
{
  boost::filesystem::ofstream output_stream(output_dir / "objectives");
  output_stream << objectives.size() << "\n";
  set_stream_precision(output_stream);
  for(auto &objective : objectives)
    {
      output_stream << objective << "\n";
    }
}
