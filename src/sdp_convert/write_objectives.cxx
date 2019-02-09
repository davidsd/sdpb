#include "write_vector.hxx"
#include "../set_stream_precision.hxx"

void write_objectives(const boost::filesystem::path &output_dir,
                      const El::BigFloat &objective_const,
                      const std::vector<El::BigFloat> &dual_objective_b)
{
  boost::filesystem::ofstream output_stream(output_dir / "objectives");
  set_stream_precision(output_stream);
  output_stream << objective_const << "\n";
  write_vector(output_stream, dual_objective_b);
}
