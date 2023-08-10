#include "diff.hxx"

#include "Float.hxx"
#include "json.hxx"

#include "sdp_convert/Dual_Constraint_Group.hxx"

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

  struct Parse_Block
  {
    Dual_Constraint_Group group;

    explicit Parse_Block(const boost::filesystem::path &block_path_no_extension)
    {
      CAPTURE(block_path_no_extension);

      auto block_path = change_extension(block_path_no_extension, ".bin");
      if(exists(block_path))
        {
          group = parse_bin(block_path);
        }
      else
        {
          block_path = change_extension(block_path_no_extension, ".json");
          REQUIRE(exists(block_path));
          group = parse_json(block_path);
        }
    }

  private:
    static Dual_Constraint_Group
    parse_json(const boost::filesystem::path &block_path)
    {
      CAPTURE(block_path);
      boost::filesystem::ifstream is(block_path);
      rapidjson::IStreamWrapper wrapper(is);
      rapidjson::Document document;
      document.ParseStream(wrapper);

      Dual_Constraint_Group result;
      result.dim = document["dim"].GetInt();
      result.num_points = document["num_points"].GetInt();
      result.bilinear_bases[0]
        = parse_Float_Matrix(document["bilinear_bases_even"]);
      result.bilinear_bases[1]
        = parse_Float_Matrix(document["bilinear_bases_odd"]);
      result.constraint_constants = parse_Float_Vector(document["c"]);
      result.constraint_matrix = parse_Float_Matrix(document["B"]);
      return result;
    }
    static Dual_Constraint_Group
    parse_bin(const boost::filesystem::path &block_path)
    {
      CAPTURE(block_path);
      boost::filesystem::ifstream is(block_path, std::ios::binary);
      boost::archive::binary_iarchive ar(is);
      Dual_Constraint_Group result;
      ar >> result;
      return result;
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
    DIFF(a.num_blocks, b.num_blocks);
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
    DIFF(a.constant, b.constant);
    DIFF(a.b, b.b);
  }

  void
  diff_block_file(const boost::filesystem::path &a_block_path_no_extension,
                  const boost::filesystem::path &b_block_path_no_extension)
  {
    INFO("diff block files");
    CAPTURE(a_block_path_no_extension);
    CAPTURE(b_block_path_no_extension);

    Parse_Block a(a_block_path_no_extension);
    Parse_Block b(b_block_path_no_extension);

    DIFF(a.group.dim, b.group.dim);
    DIFF(a.group.num_points, b.group.num_points);
    DIFF(a.group.bilinear_bases[0], b.group.bilinear_bases[0]);
    DIFF(a.group.bilinear_bases[1], b.group.bilinear_bases[1]);
    DIFF(a.group.constraint_constants, b.group.constraint_constants);
    DIFF(a.group.constraint_matrix, b.group.constraint_matrix);
  }
}

// Implementation
namespace Test_Util::REQUIRE_Equal
{
  void diff_sdp_zip(const boost::filesystem::path &a_sdp_zip,
                    const boost::filesystem::path &b_sdp_zip,
                    unsigned int input_precision, unsigned int diff_precision,
                    Test_Case_Runner runner)
  {
    INFO("diff sdp.zip files");
    CAPTURE(a_sdp_zip);
    CAPTURE(b_sdp_zip);
    Float_Binary_Precision prec(input_precision, diff_precision);
    auto a = runner.unzip_to_temp_dir(a_sdp_zip);
    auto b = runner.unzip_to_temp_dir(b_sdp_zip);
    diff_control_json(a / "control.json", b / "control.json");
    diff_objectives_json(a / "objectives.json", b / "objectives.json");
    Parse_Control_Json control(a / "control.json");
    for(int i = 0; i < control.num_blocks; ++i)
      {
        auto block_name = "block_" + std::to_string(i);
        diff_block_file(a / block_name, b / block_name);
      }
  }
}
