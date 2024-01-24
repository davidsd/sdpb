#include "pmp_read/read_json/Json_PMP_Parser.hxx"

#include <rapidjson/istreamwrapper.h>
#include <rapidjson/error/en.h>

namespace fs = std::filesystem;

PMP_File_Parse_Result
read_json(const std::filesystem::path &input_path,
          const std::function<bool(size_t matrix_index)> &should_parse_matrix)
{
  std::ifstream input_file(input_path);
  rapidjson::IStreamWrapper wrapper(input_file);
  PMP_File_Parse_Result result;
  Json_PMP_Parser parser(
    should_parse_matrix,
    [&](PMP_File_Parse_Result &&value) { result = std::move(value); });

  rapidjson::ParseResult res;
  try
    {
      rapidjson::Reader reader;
      res = reader.Parse(wrapper, parser);
    }
  catch(std::exception &e)
    {
      El::RuntimeError("Failed to parse ", input_path,
                       ": offset=", wrapper.Tell(), ": ", e.what());
    }
  if(res.IsError())
    {
      El::RuntimeError("Failed to parse ", input_path,
                       ": offset=", res.Offset(),
                       ": error: ", rapidjson::GetParseError_En(res.Code()));
    }
  return result;
}
