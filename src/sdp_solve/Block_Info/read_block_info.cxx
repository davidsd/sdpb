#include "../Block_Info.hxx"

#include <rapidjson/reader.h>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <boost/filesystem/fstream.hpp>

namespace
{
  struct Block_Info_Parser
      : public rapidjson::BaseReaderHandler<rapidjson::UTF8<>,
                                            Block_Info_Parser>
  {
    bool parsing_dim = false, parsing_num_points = false;
    int64_t dim = 0, num_points = 0;

    bool Null() { return true; }
    bool Bool(bool) { return true; }
    bool parse_integer(const int64_t &value)
    {
      if(parsing_dim)
        {
          dim = value;
          parsing_dim = false;
        }
      else if(parsing_num_points)
        {
          num_points = value;
          parsing_num_points = false;
        }
      // Quit early if finished reading dim and num_points
      return (dim == 0 || num_points == 0);
    }
    bool Int(int value) { return parse_integer(value); }
    bool Uint(unsigned value) { return parse_integer(value); }
    bool Int64(int64_t value) { return parse_integer(value); }
    bool Uint64(uint64_t value) { return parse_integer(value); }
    bool Double(double) { return true; }
    bool RawNumber(const Ch *, rapidjson::SizeType, bool) { return true; }
    bool String(const Ch *, rapidjson::SizeType, bool) { return true; }
    bool StartObject() { return true; }
    bool Key(const Ch *str, rapidjson::SizeType length, bool)
    {
      const std::string key(str, length);
      if("dim" == key)
        {
          parsing_dim = true;
        }
      else if("num_points" == key)
        {
          parsing_num_points = true;
        }
      return true;
    }
    bool EndObject(rapidjson::SizeType) { return true; }
    bool StartArray() { return true; }
    bool EndArray(rapidjson::SizeType) { return true; }
  };
}

void Block_Info::read_block_info(const boost::filesystem::path &sdp_directory)
{
  const size_t num_blocks([&]() {
    boost::filesystem::ifstream control_stream(sdp_directory / "control.json");
    rapidjson::IStreamWrapper wrapper(control_stream);
    rapidjson::Document d;
    d.ParseStream(wrapper);
    return d["num_blocks"].GetInt();
  }());

  for(size_t block(0); block != num_blocks; ++block)
    {
      boost::filesystem::path block_path(
        sdp_directory / ("block_" + std::to_string(block) + ".json"));
      boost::filesystem::ifstream block_stream(block_path);

      rapidjson::IStreamWrapper wrapper(block_stream);
      Block_Info_Parser parser;
      rapidjson::Reader reader;
      reader.Parse(wrapper, parser);
      if(parser.dim == 0 || parser.num_points == 0)
        {
          throw std::runtime_error("Unable to parse " + block_path.native());
        }
      dimensions.push_back(parser.dim);
      num_points.push_back(parser.num_points);
    }
}
