#include "Json_Function_Blocks_Parser.hxx"
#include "outer_limits/Function.hxx"
#include "pmp_read/pmp_read.hxx"

#include <filesystem>
#include <rapidjson/istreamwrapper.h>

namespace fs = std::filesystem;

void read_json(
  const fs::path &input_path, std::vector<El::BigFloat> &objectives,
  std::vector<El::BigFloat> &normalization,
  std::vector<std::vector<std::vector<std::vector<Function>>>> &functions);

void read_function_blocks(
  const fs::path &input_file, std::vector<El::BigFloat> &objectives,
  std::vector<El::BigFloat> &normalization,
  std::vector<std::vector<std::vector<std::vector<Function>>>> &functions)
{
  for(auto &json_file : collect_files_expanding_nsv(input_file))
    {
      try
        {
          ASSERT(json_file.extension() == ".json", json_file);
          std::ifstream ifs(json_file);
          rapidjson::IStreamWrapper wrapper(ifs);

          Json_Functions_Parser parser(
            [&objectives, &normalization,
             &functions](Functions_File_Parse_Result &&result) {
              objectives = std::move(result.objective);
              normalization = std::move(result.normalization);
              functions = std::move(result.functions);
            });

          rapidjson::Reader reader;
          reader.Parse(wrapper, parser);
        }
      catch(std::exception &e)
        {
          RUNTIME_ERROR("Failed to parse ", json_file, "\n", e.what());
        }
    }
}
