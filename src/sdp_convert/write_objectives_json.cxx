#include "write_vector.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <El.hpp>
#include <iostream>

void write_objectives_json(std::ostream &output_stream,
                           const El::BigFloat &objective_const,
                           const std::vector<El::BigFloat> &dual_objective_b)
{
  set_stream_precision(output_stream);
  output_stream << "{\n  \"constant\": \"" << objective_const << "\",\n";
  write_vector(output_stream, dual_objective_b, "  ", "b");
  output_stream << "\n}\n";
}
