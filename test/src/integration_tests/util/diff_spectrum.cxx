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
  struct Zero
  {
    Float zero;
    Float_Vector lambda;
  };
  struct Zeros
  {
    std::vector<Zero> zeros;
    Float error;
  };

  struct Parse_Spectrum_Json : boost::noncopyable
  {
    std::vector<Zeros> zeros_array;

    std::vector<std::pair<std::string, std::string>> options;
    explicit Parse_Spectrum_Json(const boost::filesystem::path &path,
                                 unsigned int binary_precision)
    {
      CAPTURE(path);
      REQUIRE(exists(path));
      Float_Binary_Precision _(binary_precision);
      boost::filesystem::ifstream is(path);
      rapidjson::IStreamWrapper wrapper(is);
      rapidjson::Document document;
      document.ParseStream(wrapper);

      for(const auto &item : document.GetArray())
        {
          Zeros zeros;
          for(const auto &z : item["zeros"].GetArray())
            {
              Zero zero;
              zero.zero = Json::parse_Float(z["zero"]);
              zero.lambda = Json::parse_Float_Vector(z["lambda"]);
              zeros.zeros.emplace_back(zero);
            }
          zeros.error = Json::parse_Float(item["error"]);
          zeros_array.emplace_back(zeros);
        }
    }
  };

  using Test_Util::REQUIRE_Equal::diff;

  void diff(const Zero &a, const Zero &b)
  {
    diff(a.zero, b.zero);
    diff(a.lambda, b.lambda);
  }
  void diff(const Zeros &a, const Zeros &b)
  {
    diff(a.error, b.error);
    diff(a.zeros, b.zeros);
  }
}

// Implementation
namespace Test_Util::REQUIRE_Equal
{
  void diff_spectrum(const boost::filesystem::path &a_json,
                     const boost::filesystem::path &b_json,
                     unsigned int binary_precision)
  {
    INFO("diff spectrum output");
    CAPTURE(a_json);
    CAPTURE(b_json);
    Parse_Spectrum_Json a(a_json, binary_precision);
    Parse_Spectrum_Json b(b_json, binary_precision);
    diff(a.zeros_array, b.zeros_array);
  }
}
