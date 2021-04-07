#include "../Block_Info.hxx"

#include <archive_reader.hpp>
#include <archive_exception.hpp>

#include <rapidjson/reader.h>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>

#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string/predicate.hpp>

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

  void
  parse_block(const size_t &block_number, std::istream &block_stream,
              std::vector<size_t> &dimensions, std::vector<size_t> &num_points)
  {
    rapidjson::IStreamWrapper wrapper(block_stream);
    Block_Info_Parser parser;
    rapidjson::Reader reader;
    reader.Parse(wrapper, parser);
    if(parser.dim == 0 || parser.num_points == 0)
      {
        throw std::runtime_error("Unable to parse block: "
                                 + std::to_string(block_number));
      }
    dimensions.at(block_number) = parser.dim;
    num_points.at(block_number) = parser.num_points;
  }

  size_t parse_num_blocks(std::istream &control_stream)
  {
    rapidjson::IStreamWrapper wrapper(control_stream);
    rapidjson::Document document;
    document.ParseStream(wrapper);
    return document["num_blocks"].GetInt();
  }
}

void Block_Info::read_block_info(const boost::filesystem::path &sdp_path)
{
  if(boost::filesystem::is_regular_file(sdp_path))
    {
      const size_t num_blocks([&]() {
        boost::filesystem::ifstream fs(sdp_path);
        ns_archive::reader reader = ns_archive::reader::make_reader<
          ns_archive::ns_reader::format::_ALL,
          ns_archive::ns_reader::filter::_ALL>(fs, 10240);

        for(auto entry : reader)
          {
            if(entry->get_header_value_pathname() == "control.json")
              {
                return parse_num_blocks(entry->get_stream());
              }
          }
        throw std::runtime_error(
          "Unable to find control.json in sdp input file");
      }());

      dimensions.resize(num_blocks);
      num_points.resize(num_blocks);

      boost::filesystem::ifstream fs(sdp_path);
      ns_archive::reader reader
        = ns_archive::reader::make_reader<ns_archive::ns_reader::format::_ALL,
                                          ns_archive::ns_reader::filter::_ALL>(
          fs, 10240);

      const std::string prefix("block_");
      for(auto entry : reader)
        {
          const std::string pathname(entry->get_header_value_pathname());
          if(boost::algorithm::starts_with(pathname, prefix))
            {
              const size_t block_number(
                std::stoll(pathname.substr(prefix.size())));
              if(block_number >= num_blocks)
                {
                  throw std::runtime_error(
                    "Invalid block number for entry '" + pathname + "' in '"
                    + sdp_path.string()
                    + "'. The block number must be between 0 and "
                    + std::to_string(num_blocks - 1) + ".");
                }

              parse_block(block_number, entry->get_stream(), dimensions,
                          num_points);
            }
        }
      for(size_t block_number(0); block_number != dimensions.size();
          ++block_number)
        {
          if(dimensions.at(block_number) == 0)
            {
              throw std::runtime_error(
                "Missing block " + std::to_string(block_number)
                + " from sdp path: '" + sdp_path.string() + "'.");
            }
        }
    }
  else
    {
      const size_t num_blocks([&]() {
        boost::filesystem::ifstream control_stream(sdp_path / "control.json");
        return parse_num_blocks(control_stream);
      }());

      dimensions.resize(num_blocks);
      num_points.resize(num_blocks);
      for(size_t block(0); block != num_blocks; ++block)
        {
          boost::filesystem::path block_path(
            sdp_path / ("block_" + std::to_string(block) + ".json"));
          boost::filesystem::ifstream block_stream(block_path);
          parse_block(block, block_stream, dimensions, num_points);
        }
    }
}
