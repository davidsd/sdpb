#include "../Block_Info.hxx"
#include "../Archive_Reader.hxx"
#include "../../sdp_convert.hxx"

#include <rapidjson/reader.h>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>

#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string/predicate.hpp>

namespace
{
  void parse_block_info_json(const size_t &block_index,
                             std::istream &block_info_stream,
                             std::vector<size_t> &dimensions,
                             std::vector<size_t> &num_points)
  {
    rapidjson::IStreamWrapper wrapper(block_info_stream);
    rapidjson::Document document;
    document.ParseStream(wrapper);
    size_t dim_value = document["dim"].GetInt64();
    size_t num_points_value = document["num_points"].GetInt64();
    if(dim_value == 0 || num_points_value == 0)
      {
        El::RuntimeError("Unable to parse block_", block_index,
                         ": dim=", dim_value,
                         ", num_points=", num_points_value);
      }

    dimensions.at(block_index) = dim_value;
    num_points.at(block_index) = num_points_value;
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
        Archive_Reader reader(sdp_path);
        while(reader.next_entry())
          {
            using namespace std::string_literals;
            if("control.json"s == archive_entry_pathname(reader.entry_ptr))
              {
                std::istream stream(&reader);
                return parse_num_blocks(stream);
              }
          }
        El::RuntimeError("Unable to find control.json in sdp input file");
      }());

      dimensions.resize(num_blocks);
      num_points.resize(num_blocks);

      const std::string prefix("block_info_");
      Archive_Reader reader(sdp_path);
      while(reader.next_entry())
        {
          const std::string pathname(archive_entry_pathname(reader.entry_ptr));
          if(boost::algorithm::starts_with(pathname, prefix))
            {
              const size_t block_index(
                std::stoll(pathname.substr(prefix.size())));
              if(block_index >= num_blocks)
                {
                  El::RuntimeError("Invalid block number for entry '",
                                   pathname, "' in '", sdp_path,
                                   ". The block number must be between 0 and ",
                                   num_blocks - 1, ".");
                }
              std::istream stream(&reader);
              parse_block_info_json(block_index, stream, dimensions,
                                    num_points);
            }
        }
      for(size_t block_index(0); block_index != num_blocks; ++block_index)
        {
          if(dimensions.at(block_index) == 0
             || num_points.at(block_index) == 0)
            {
              El::RuntimeError("Missing block ", block_index,
                               " from sdp path: ", sdp_path);
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
          boost::filesystem::path block_info_path(
            sdp_path / ("block_info_" + std::to_string(block) + ".json"));
          boost::filesystem::ifstream block_stream(block_info_path);
          parse_block_info_json(block, block_stream, dimensions, num_points);
        }
    }
}
