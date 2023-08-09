#include "diff.hxx"

#include "Float.hxx"
#include "json.hxx"

#include <catch2/catch_amalgamated.hpp>

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>

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
    explicit Parse_Outer_Limits_Json(const boost::filesystem::path &path)
    {
      CAPTURE(path);
      REQUIRE(exists(path));
      boost::filesystem::ifstream is(path);
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
  diff_outer_limits(const boost::filesystem::path &a_json,
                    const boost::filesystem::path &b_json,
                    unsigned int input_precision, unsigned int diff_precision)
  {
    INFO("diff outer_limits output");
    CAPTURE(a_json);
    CAPTURE(b_json);
    Float_Binary_Precision prec(input_precision, diff_precision);
    Parse_Outer_Limits_Json a(a_json);
    Parse_Outer_Limits_Json b(b_json);
    diff(a.optimal, b.optimal);
    diff(a.y, b.y);
    diff(a.options, b.options);
  }
}
