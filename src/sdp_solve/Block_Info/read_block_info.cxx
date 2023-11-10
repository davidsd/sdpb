#include "../Block_Info.hxx"
#include "../Archive_Reader.hxx"
#include "../../sdp_convert.hxx"

#include <rapidjson/reader.h>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>

#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <iostream> 

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

    template <typename T> bool parse_integer(const T &value)
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
      else
        {
          throw std::runtime_error("Invalid input file.  Found the integer '"
                                   + std::to_string(value)
                                   + "' outside of 'dim' or 'num_points'.");
        }
      return true;
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

  Block_File_Format get_block_format(const boost::filesystem::path &block_path)
  {
    auto extension = block_path.extension();
    if(extension == ".bin")
      return Block_File_Format::bin;
    else if(extension == ".json")
      return Block_File_Format::json;
    El::RuntimeError("Unknown block file extension: ", block_path);
  }

  void parse_block(const size_t &block_index, std::istream &block_stream,
                   Block_File_Format format, std::vector<size_t> &dimensions,
                   std::vector<size_t> &num_points)
  {
    size_t dim_value;
    size_t num_points_value;
    if(format == bin)
      {
        // TODO reads the whole file and allocates unnecessary RAM
        // We need to read all bytes anyway,
        // otherwise libarchive fails to reasd next entry
        // TODO store dim and num_points in a separate file
        boost::archive::binary_iarchive ar(block_stream);
        Dual_Constraint_Group group;
        ar >> group;
        dim_value = group.dim;
        num_points_value = group.num_points;
        if(group.block_index != block_index)
          {
            El::RuntimeError("Read wrong block_index=", group.block_index,
                             " from block_", block_index);
          }
      }
    else if(format == json)
      {
        rapidjson::IStreamWrapper wrapper(block_stream);
        Block_Info_Parser parser;
        rapidjson::Reader reader;
        reader.Parse(wrapper, parser);
        dim_value = parser.dim;
        num_points_value = parser.num_points;
      }
    else
      El::RuntimeError("Unknown SDP format: ", format);

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
        throw std::runtime_error(
          "Unable to find control.json in sdp input file");
      }());

      dimensions.resize(num_blocks);
      num_points.resize(num_blocks);

      const std::string prefix("block_");
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
                  throw std::runtime_error(
                    "Invalid block number for entry '" + pathname + "' in '"
                    + sdp_path.string()
                    + "'. The block number must be between 0 and "
                    + std::to_string(num_blocks - 1) + ".");
                }
              Block_File_Format format = get_block_format(pathname);
              std::istream stream(&reader);
              parse_block(block_index, stream, format, dimensions, num_points);
            }
        }
      for(size_t block_index(0); block_index != dimensions.size();
          ++block_index)
        {
          if(dimensions.at(block_index) == 0)
            {
              throw std::runtime_error(
                "Missing block " + std::to_string(block_index)
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
            sdp_path / ("block_" + std::to_string(block) + ".bin"));
          if(!exists(block_path))
            block_path = change_extension(block_path, ".json");
          Block_File_Format format = get_block_format(block_path);
          boost::filesystem::ifstream block_stream(block_path,
                                                   std::ios::binary);
          parse_block(block, block_stream, format, dimensions, num_points);
        }
    }
}
