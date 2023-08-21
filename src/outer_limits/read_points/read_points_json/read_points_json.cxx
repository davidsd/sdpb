#include "Points_Parser.hxx"

#include <rapidjson/istreamwrapper.h>
#include <filesystem>

namespace fs = std::filesystem;

void read_points_json(const fs::path &input_path,
                      std::vector<std::vector<El::BigFloat>> &points)
{
  std::ifstream input_file(input_path);
  rapidjson::IStreamWrapper wrapper(input_file);
  Points_Parser parser;
  rapidjson::Reader reader;
  reader.Parse(wrapper, parser);

  for(auto &block_points: parser.points_state.value)
    {
      points.push_back(block_points);
    }
}
