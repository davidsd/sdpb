#include "diff.hxx"

#include "Float.hxx"
#include "json.hxx"

#include <catch2/catch_amalgamated.hpp>

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>

using namespace Test_Util::Json;

// Parser classes
// TODO replace them with faster parsers from SDPB
// These parsers use DOM, which is simpler but slower in comparison to SAX
namespace
{
  struct Parse_Control_Json : boost::noncopyable
  {
    int num_blocks;
    std::string command;
    explicit Parse_Control_Json(const boost::filesystem::path &path)
    {
      CAPTURE(path);
      REQUIRE(exists(path));
      boost::filesystem::ifstream is(path);
      rapidjson::IStreamWrapper wrapper(is);
      rapidjson::Document document;
      document.ParseStream(wrapper);

      num_blocks = document["num_blocks"].GetInt();
      command = document["command"].GetString();
    }
  };

  struct Parse_Objectives_Json : boost::noncopyable
  {
    Float constant;
    std::vector<Float> b;
    explicit Parse_Objectives_Json(const boost::filesystem::path &path)
    {
      CAPTURE(path);
      REQUIRE(exists(path));
      boost::filesystem::ifstream is(path);
      rapidjson::IStreamWrapper wrapper(is);
      rapidjson::Document document;
      document.ParseStream(wrapper);

      constant = parse_Float(document["constant"]);
      b = parse_Float_Vector(document["b"]);
    }
  };

  struct Parse_Block_Json : boost::noncopyable
  {
    int dim;
    int num_points;
    Float_Matrix bilinear_bases_even;
    Float_Matrix bilinear_bases_odd;
    Float_Vector c;
    Float_Matrix B;

    explicit Parse_Block_Json(const boost::filesystem::path &path)
    {
      CAPTURE(path);
      REQUIRE(exists(path));
      boost::filesystem::ifstream is(path);
      rapidjson::IStreamWrapper wrapper(is);
      rapidjson::Document document;
      document.ParseStream(wrapper);

      dim = document["dim"].GetInt();
      num_points = document["num_points"].GetInt();
      bilinear_bases_even
        = parse_Float_Matrix(document["bilinear_bases_even"]);
      bilinear_bases_odd = parse_Float_Matrix(document["bilinear_bases_odd"]);
      c = parse_Float_Vector(document["c"]);
      B = parse_Float_Matrix(document["B"]);
    }
  };
}

// Helper functions
namespace
{
  using Test_Util::REQUIRE_Equal::diff;

  void diff_control_json(const boost::filesystem::path &a_control_json,
                         const boost::filesystem::path &b_control_json)
  {
    CAPTURE(a_control_json);
    CAPTURE(b_control_json);
    Parse_Control_Json a(a_control_json);
    Parse_Control_Json b(b_control_json);
    diff(a.num_blocks, b.num_blocks);
    // ignore "command", since it's unimportant:
    // diff(a.command, b.command);
  }

  void diff_objectives_json(const boost::filesystem::path &a_objectives_json,
                            const boost::filesystem::path &b_objectives_json)
  {
    CAPTURE(a_objectives_json);
    CAPTURE(b_objectives_json);
    Parse_Objectives_Json a(a_objectives_json);
    Parse_Objectives_Json b(b_objectives_json);
    diff(a.constant, b.constant);
    diff(a.b, b.b);
  }

  void diff_block_json(const boost::filesystem::path &a_block_json,
                       const boost::filesystem::path &b_block_json)

  {
    INFO("diff_block_json");
    CAPTURE(a_block_json);
    CAPTURE(b_block_json);
    Parse_Block_Json a_block(a_block_json);
    Parse_Block_Json b_block(b_block_json);

    INFO("diff dim");
    diff(a_block.dim, b_block.dim);
    INFO("diff num_points");
    diff(a_block.num_points, b_block.num_points);
    INFO("diff bilinear_bases_even");
    diff(a_block.bilinear_bases_even, b_block.bilinear_bases_even);
    INFO("diff bilinear_bases_odd");
    diff(a_block.bilinear_bases_odd, b_block.bilinear_bases_odd);
    INFO("diff constraint vector c");
    diff(a_block.c, b_block.c);
    INFO("diff constraint matrix B");
    diff(a_block.B, b_block.B);
  }
}

// Implementation
namespace Test_Util::REQUIRE_Equal
{
  void diff_sdp_zip(const boost::filesystem::path &a_sdp_zip,
                    const boost::filesystem::path &b_sdp_zip,
                    unsigned int binary_precision, Test_Case_Runner runner)
  {
    INFO("diff sdp.zip files");
    CAPTURE(a_sdp_zip);
    CAPTURE(b_sdp_zip);
    Float_Binary_Precision prec(binary_precision);
    auto a = runner.unzip_to_temp_dir(a_sdp_zip);
    auto b = runner.unzip_to_temp_dir(b_sdp_zip);
    diff_control_json(a / "control.json", b / "control.json");
    diff_objectives_json(a / "objectives.json", b / "objectives.json");
    Parse_Control_Json control(a / "control.json");
    for(int i = 0; i < control.num_blocks; ++i)
      {
        auto block_name = "block_" + std::to_string(i) + ".json";
        diff_block_json(a / block_name, b / block_name);
      }
  }
}
