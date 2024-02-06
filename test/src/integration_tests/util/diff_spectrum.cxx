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
  struct Zero
  {
    Float zero;
    Float_Vector lambda;
  };
  struct Zeros
  {
    std::vector<Zero> zeros;
    Float error;
    std::string block_path;
  };

  struct Parse_Spectrum_Json : boost::noncopyable
  {
    std::vector<Zeros> zeros_array;

    std::vector<std::pair<std::string, std::string>> options;
    explicit Parse_Spectrum_Json(const fs::path &path)
    {
      CAPTURE(path);
      REQUIRE(exists(path));
      std::ifstream is(path);
      rapidjson::IStreamWrapper wrapper(is);
      rapidjson::Document document;
      document.ParseStream(wrapper);

      for(const auto &item : document.GetArray())
        {
          Zeros zeros;
          zeros.block_path = item["block_path"].GetString();
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
    DIFF(a.zero, b.zero);
    DIFF(a.lambda, b.lambda);
  }
  void diff(const Zeros &a, const Zeros &b)
  {
    DIFF(a.block_path, b.block_path);
    DIFF(a.error, b.error);
    DIFF(a.zeros, b.zeros);
  }
}

// Implementation
namespace Test_Util::REQUIRE_Equal
{
  void diff_spectrum(const fs::path &a_json, const fs::path &b_json,
                     unsigned int input_precision, unsigned int diff_precision)
  {
    INFO("diff spectrum output");
    CAPTURE(a_json);
    CAPTURE(b_json);
    Float_Binary_Precision prec(input_precision, diff_precision);
    Parse_Spectrum_Json a(a_json);
    Parse_Spectrum_Json b(b_json);
    DIFF(a.zeros_array, b.zeros_array);
  }
}
