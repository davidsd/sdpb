#include "diff.hxx"

#include "Float.hxx"
#include "json.hxx"

#include <catch2/catch_amalgamated.hpp>

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>

namespace fs = std::filesystem;
using namespace Test_Util;

// Parser classes
// TODO replace them with faster parsers from SDPB
// These parsers use DOM, which is simpler but slower in comparison to SAX
namespace
{
  struct Parse_Outer_Limits_Json : boost::noncopyable
  {
    Float optimal;
    Float_Vector y;
    std::vector<std::pair<std::string, std::string>> options;
    explicit Parse_Outer_Limits_Json(const fs::path &path)
    {
      CAPTURE(path);
      REQUIRE(exists(path));
      std::ifstream is(path);
      rapidjson::IStreamWrapper wrapper(is);
      rapidjson::Document document;
      document.ParseStream(wrapper);

      optimal = Json::parse_Float(document["optimal"]);
      y = Json::parse_Float_Vector(document["y"]);
      for(const auto &option : document["options"].GetObject())
        {
          options.emplace_back(option.name.GetString(),
                               option.value.GetString());
        }
    }
  };
}

// Implementation
namespace Test_Util::REQUIRE_Equal
{
  void
  diff_outer_limits(const fs::path &a_json, const fs::path &b_json,
                    unsigned int input_precision, unsigned int diff_precision)
  {
    INFO("diff outer_limits output");
    CAPTURE(a_json);
    CAPTURE(b_json);
    Float_Binary_Precision prec(input_precision, diff_precision);
    Parse_Outer_Limits_Json a(a_json);
    Parse_Outer_Limits_Json b(b_json);
    DIFF(a.optimal, b.optimal);
    DIFF(a.y, b.y);
    DIFF(a.options.size(), b.options.size());
    for(size_t i = 0; i < a.options.size(); ++i)
      {
        const auto &[a_key, a_value] = a.options.at(i);
        const auto &[b_key, b_value] = b.options.at(i);
        DIFF(a_key, b_key);
        // Do not compare output paths
        if(a_key == "initialCheckpointDir" || a_key == "checkpointDir"
           || a_key == "out")
          continue;
        CAPTURE(a_key);
        DIFF(a_value, b_value);
      }
  }
}
