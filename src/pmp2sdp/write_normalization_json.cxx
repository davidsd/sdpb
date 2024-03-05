#include "write_vector.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <El.hpp>
#include <iostream>

void write_normalization_json(std::ostream &output_stream,
                           const std::vector<El::BigFloat> &normalization)
{
  set_stream_precision(output_stream);
  output_stream << "{\n";
  write_vector(output_stream, normalization, "  ", "normalization");
  output_stream << "\n}\n";
}
