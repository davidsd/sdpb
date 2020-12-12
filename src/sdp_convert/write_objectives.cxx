#include "write_vector.hxx"
#include "../set_stream_precision.hxx"

void write_objectives(const boost::filesystem::path &output_dir,
                      const El::BigFloat &objective_const,
                      const std::vector<El::BigFloat> &dual_objective_b)
{
  const boost::filesystem::path output_path(output_dir / "objectives.json");
  boost::filesystem::ofstream output_stream(output_path);
  set_stream_precision(output_stream);
  output_stream << "{\n  \"c\": \"" << objective_const << "\",\n  ";
  write_vector(output_stream, dual_objective_b, "b");
  output_stream << "\n}\n";
  if(!output_stream.good())
    {
      throw std::runtime_error("Error when writing to: "
                               + output_path.string());
    }
}
