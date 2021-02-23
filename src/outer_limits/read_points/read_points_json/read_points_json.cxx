#include "Points_Parser.hxx"

#include <rapidjson/istreamwrapper.h>
#include <boost/filesystem/fstream.hpp>

void read_points_json(const boost::filesystem::path &input_path,
                      std::vector<std::vector<El::BigFloat>> &points)
{
  boost::filesystem::ifstream input_file(input_path);
  rapidjson::IStreamWrapper wrapper(input_file);
  Points_Parser parser;
  rapidjson::Reader reader;
  reader.Parse(wrapper, parser);

  for(auto &block_points: parser.points_state.value)
    {
      points.push_back(block_points);
    }
}
