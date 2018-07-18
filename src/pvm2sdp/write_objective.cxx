#include "../set_stream_precision.hxx"
#include "write_vector.hxx"

void write_objective(const boost::filesystem::path &output_dir,
                     const std::vector<El::BigFloat> &objectives)
{
  boost::filesystem::ofstream output_stream(output_dir / "objectives");
  set_stream_precision(output_stream);
  write_vector(output_stream,objectives);
}
