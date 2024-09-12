#include "pmp/PMP_Info.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/json/Abstract_Json_Object_Parser.hxx"
#include "sdpb_util/json/Json_Float_Parser.hxx"
#include "sdpb_util/json/Json_Vector_Parser.hxx"
#include "sdpb_util/json/Json_Vector_Parser_With_Skip.hxx"
#include "sdpb_util/json/parse_json.hxx"

#include <El.hpp>

#include <filesystem>
#include <set>

using Json_Block_Vector_Parser
  = Json_Vector_Parser_With_Skip<Json_Vector_Parser<Json_BigFloat_Parser>>;

class Json_c_minus_By_Parser final
    : public Abstract_Json_Object_Parser<std::vector<El::Matrix<El::BigFloat>>>
{
private:
  // Each block in result is a vector represented as matrix with width = 1
  std::vector<El::Matrix<El::BigFloat>> result;
  Json_Block_Vector_Parser inner_parser;

public:
  Json_c_minus_By_Parser(
    const std::function<bool(size_t matrix_index)> &should_skip_block,
    const std::function<void(value_type &&result)> &on_parsed)
      : Abstract_Json_Object_Parser(false, on_parsed, [] {}),
        inner_parser(
          false,
          [this](Vector_Parse_Result_With_Skip<std::vector<El::BigFloat>>
                   &&parse_result) {
            for(const auto &vec : parse_result.parsed_elements)
              {
                const int height = vec.size();
                const int width = 1;
                auto &m = this->result.emplace_back(height, width);
                for(int i = 0; i < height; ++i)
                  m.Set(i, 0, vec.at(i));
              }
          },
          [] {}, should_skip_block)
  {}
  value_type get_result() override { return result; }

protected:
  Abstract_Json_Reader_Handler &element_parser(const std::string &key) override
  {
    if(key == "c_minus_By")
      return inner_parser;
    // TODO print warning instead?
    RUNTIME_ERROR("Unknown key: ", key);
  }

public:
  void reset_element_parsers(bool skip) override { inner_parser.reset(skip); }
  void clear_result() override { result.clear(); }
};

std::vector<El::Matrix<El::BigFloat>>
read_c_minus_By(const std::filesystem::path &input_path,
                const PMP_Info &pmp_info)
{
  std::vector<El::Matrix<El::BigFloat>> c_minus_By_blocks;
  {
    std::set<size_t> block_indices;
    for(const auto &block : pmp_info.blocks)
      block_indices.emplace(block.block_index);

    const auto should_skip_block
      = [&block_indices](const size_t global_index) {
          return block_indices.find(global_index) == block_indices.end();
        };
    const auto on_parsed = [&c_minus_By_blocks](auto &&result) {
      c_minus_By_blocks = std::forward<decltype(result)>(result);
    };
    Json_c_minus_By_Parser parser(should_skip_block, on_parsed);
    parse_json(input_path, parser);
    ASSERT_EQUAL(c_minus_By_blocks.size(), block_indices.size());
    return c_minus_By_blocks;
  }
}
