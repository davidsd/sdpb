#include "pmp_read/read_json/Json_PMP_Parser.hxx"
#include "sdpb_util/json/parse_json.hxx"

namespace fs = std::filesystem;

PMP_File_Parse_Result
read_json(const std::filesystem::path &input_path, bool should_parse_objective,
          bool should_parse_normalization,
          const std::function<bool(size_t matrix_index)> &should_parse_matrix)
{
  PMP_File_Parse_Result result;
  Json_PMP_Parser parser(
    should_parse_objective, should_parse_normalization, should_parse_matrix,
    [&](PMP_File_Parse_Result &&value) { result = std::move(value); });

  parse_json(input_path, parser);
  return result;
}
