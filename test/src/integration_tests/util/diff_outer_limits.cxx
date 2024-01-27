#include "diff.hxx"

#include "Float.hxx"
#include "json.hxx"
#include "outer_limits/Function.hxx"

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
  struct Parse_Functions_Json : boost::noncopyable
  {
    struct Function
    {
      Float max_delta;
      Float epsilon_value;
      Float infinity_value;
      Float_Vector chebyshev_values;
    };

    Float_Vector objective;
    Float_Vector normalization;
    std::vector<std::vector<std::vector<std::vector<Function>>>> functions;
    explicit Parse_Functions_Json(const fs::path &path)
    {
      CAPTURE(path);
      REQUIRE(exists(path));
      std::ifstream is(path);
      rapidjson::IStreamWrapper wrapper(is);
      rapidjson::Document document;
      document.ParseStream(wrapper);

      using Json::parse_Float;
      using Json::parse_Float_Vector;
      using Json::parse_vector;

      objective = parse_Float_Vector(document["objective"]);
      normalization = parse_Float_Vector(document["normalization"]);

      auto parse_Function = [](const Json::Json_Value &json_value) {
        Function result;
        result.max_delta = Json::parse_Float(json_value["max_delta"]);
        result.epsilon_value = Json::parse_Float(json_value["epsilon_value"]);
        result.infinity_value = parse_Float(json_value["infinity_value"]);
        result.chebyshev_values
          = parse_Float_Vector(json_value["chebyshev_values"]);
        return result;
      };

      auto parse_Function_vector = [&](const Json::Json_Value &json_value) {
        return parse_vector<Function>(json_value, parse_Function);
      };
      auto parse_Function_vector_2 = [&](const Json::Json_Value &json_value) {
        return parse_vector<std::vector<Function>>(json_value,
                                                   parse_Function_vector);
      };
      auto parse_Function_vector_3 = [&](const Json::Json_Value &json_value) {
        return parse_vector<std::vector<std::vector<Function>>>(
          json_value, parse_Function_vector_2);
      };

      functions
        = Json::parse_vector<std::vector<std::vector<std::vector<Function>>>>(
          document["functions"], parse_Function_vector_3);
    }
  };

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
  template <>
  void diff(const Parse_Functions_Json::Function &a,
            const Parse_Functions_Json::Function &b)
  {
    INFO("diff Function");
    DIFF(a.max_delta, b.max_delta);
    DIFF(a.epsilon_value, b.epsilon_value);
    DIFF(a.infinity_value, b.infinity_value);
    DIFF(a.chebyshev_values, b.chebyshev_values);
  }

  void diff_functions_json(const fs::path &a_json, const fs::path &b_json,
                           unsigned int input_precision,
                           unsigned int diff_precision)
  {
    INFO("diff outer_limits output");
    CAPTURE(a_json);
    CAPTURE(b_json);
    Float_Binary_Precision prec(input_precision, diff_precision);
    Parse_Functions_Json a(a_json);
    Parse_Functions_Json b(b_json);
    DIFF(a.objective, b.objective);
    DIFF(a.normalization, b.normalization);
    DIFF(a.functions, b.functions);
  }

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
        if(a_key == "functions" || a_key == "initialCheckpointDir"
           || a_key == "checkpointDir" || a_key == "out")
          continue;
        CAPTURE(a_key);
        DIFF(a_value, b_value);
      }
  }
}
