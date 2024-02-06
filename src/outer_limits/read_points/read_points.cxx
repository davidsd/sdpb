#include "Json_Points_Parser.hxx"
#include "pmp_read/pmp_read.hxx"

#include <rapidjson/istreamwrapper.h>

#include <filesystem>

namespace fs = std::filesystem;

void read_points(const fs::path &input_path,
                 std::vector<std::vector<El::BigFloat>> &points)
{
  for(auto &json_file : collect_files_expanding_nsv(input_path))
    {
      try
        {
          ASSERT(json_file.extension() == ".json", json_file);

          std::ifstream ifs(json_file);
          rapidjson::IStreamWrapper wrapper(ifs);

          Json_Points_Parser parser(
            [&points](std::vector<std::vector<El::BigFloat>> &&result) {
              ASSERT(!result.empty());
              for(auto &block_points : result)
                {
                  points.push_back(std::move(block_points));
                }
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
