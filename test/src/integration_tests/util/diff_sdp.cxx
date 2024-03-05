#include "diff.hxx"

#include "Float.hxx"
#include "json.hxx"
#include "sdpb_util/boost_serialization.hxx"

#include <catch2/catch_amalgamated.hpp>

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

namespace fs = std::filesystem;
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
    explicit Parse_Control_Json(const fs::path &path)
    {
      CAPTURE(path);
      REQUIRE(exists(path));
      std::ifstream is(path);
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
    explicit Parse_Objectives_Json(const fs::path &path)
    {
      CAPTURE(path);
      REQUIRE(exists(path));
      std::ifstream is(path);
      rapidjson::IStreamWrapper wrapper(is);
      rapidjson::Document document;
      document.ParseStream(wrapper);

      constant = parse_Float(document["constant"]);
      b = parse_Float_Vector(document["b"]);
    }
  };

  struct Parse_Normalization_Json : boost::noncopyable
  {
    std::vector<Float> normalization;
    explicit Parse_Normalization_Json(const fs::path &path)
    {
      CAPTURE(path);
      REQUIRE(exists(path));
      std::ifstream is(path);
      rapidjson::IStreamWrapper wrapper(is);
      rapidjson::Document document;
      document.ParseStream(wrapper);

      normalization = parse_Float_Vector(document["normalization"]);
    }
  };

  struct Parse_Block_Info_Json : boost::noncopyable
  {
    size_t dim;
    size_t num_points;
    explicit Parse_Block_Info_Json(const fs::path &path)
    {
      CAPTURE(path);
      REQUIRE(exists(path));
      std::ifstream is(path);
      rapidjson::IStreamWrapper wrapper(is);
      rapidjson::Document document;
      document.ParseStream(wrapper);

      dim = document["dim"].GetInt();
      num_points = document["num_points"].GetInt();
    }
  };

  struct Parse_Block_Data
  {
    Float_Matrix bilinear_bases_even;
    Float_Matrix bilinear_bases_odd;
    Float_Vector constraint_constants;
    Float_Matrix constraint_matrix;
    explicit Parse_Block_Data(const fs::path &block_path_no_extension)
    {
      CAPTURE(block_path_no_extension);

      auto block_path = block_path_no_extension;
      block_path.replace_extension(".bin");
      if(exists(block_path))
        {
          parse_bin(block_path);
        }
      else
        {
          block_path.replace_extension(".json");
          REQUIRE(exists(block_path));
          parse_json(block_path);
        }
    }

  private:
    void parse_json(const fs::path &block_path)
    {
      CAPTURE(block_path);
      std::ifstream is(block_path);
      rapidjson::IStreamWrapper wrapper(is);
      rapidjson::Document document;
      document.ParseStream(wrapper);

      bilinear_bases_even
        = parse_Float_Matrix(document["bilinear_bases_even"]);
      bilinear_bases_odd = parse_Float_Matrix(document["bilinear_bases_odd"]);
      constraint_constants = parse_Float_Vector(document["c"]);
      constraint_matrix = parse_Float_Matrix(document["B"]);
    }
    void parse_bin(const fs::path &block_path)
    {
      CAPTURE(block_path);
      std::ifstream is(block_path, std::ios::binary);
      boost::archive::binary_iarchive ar(is);
      mp_bitcnt_t precision;
      ar >> precision;
      REQUIRE(precision == El::gmp::Precision());
      ar >> constraint_matrix;
      ar >> constraint_constants;
      ar >> bilinear_bases_even;
      ar >> bilinear_bases_odd;
    }
  };
}

// Helper functions
namespace
{
  using Test_Util::REQUIRE_Equal::diff;

  void diff_control_json(const fs::path &a_control_json,
                         const fs::path &b_control_json)
  {
    CAPTURE(a_control_json);
    CAPTURE(b_control_json);
    Parse_Control_Json a(a_control_json);
    Parse_Control_Json b(b_control_json);
    DIFF(a.num_blocks, b.num_blocks);
    // ignore "command", since it's unimportant:
    // diff(a.command, b.command);
  }

  void diff_objectives_json(const fs::path &a_objectives_json,
                            const fs::path &b_objectives_json)
  {
    CAPTURE(a_objectives_json);
    CAPTURE(b_objectives_json);
    Parse_Objectives_Json a(a_objectives_json);
    Parse_Objectives_Json b(b_objectives_json);
    DIFF(a.constant, b.constant);
    DIFF(a.b, b.b);
  }

  void diff_normalization_json(const fs::path &a_normalization_json,
                               const fs::path &b_normalization_json)
  {
    CAPTURE(a_normalization_json);
    CAPTURE(b_normalization_json);
    Parse_Normalization_Json a(a_normalization_json);
    Parse_Normalization_Json b(b_normalization_json);
    DIFF(a.normalization, b.normalization);
  }

  void diff_block_info_json(const fs::path &a_block_info_path,
                            const fs::path &b_block_info_path)
  {
    INFO("diff block_info_XXX.json files");
    CAPTURE(a_block_info_path);
    CAPTURE(b_block_info_path);

    Parse_Block_Info_Json a(a_block_info_path);
    Parse_Block_Info_Json b(b_block_info_path);

    DIFF(a.dim, b.dim);
    DIFF(a.num_points, b.num_points);
  }
  void diff_block_data_file(const fs::path &a_block_path_no_extension,
                            const fs::path &b_block_path_no_extension)
  {
    INFO("diff block_data_XXX files");
    CAPTURE(a_block_path_no_extension);
    CAPTURE(b_block_path_no_extension);

    Parse_Block_Data a(a_block_path_no_extension);
    Parse_Block_Data b(b_block_path_no_extension);

    DIFF(a.bilinear_bases_even, b.bilinear_bases_even);
    DIFF(a.bilinear_bases_odd, b.bilinear_bases_odd);
    DIFF(a.constraint_constants, b.constraint_constants);
    DIFF(a.constraint_matrix, b.constraint_matrix);
  }
}

// Implementation
namespace Test_Util::REQUIRE_Equal
{
  void diff_sdp(const fs::path &a_sdp, const fs::path &b_sdp,
                unsigned int input_precision, unsigned int diff_precision,
                Test_Case_Runner runner, bool check_normalization)
  {
    INFO("diff sdp directories/zip files");
    REQUIRE(a_sdp != b_sdp);
    CAPTURE(a_sdp);
    CAPTURE(b_sdp);
    Float_Binary_Precision prec(input_precision, diff_precision);
    auto a = fs::is_directory(a_sdp) ? a_sdp : runner.unzip_to_temp_dir(a_sdp);
    auto b = fs::is_directory(b_sdp) ? b_sdp : runner.unzip_to_temp_dir(b_sdp);
    diff_control_json(a / "control.json", b / "control.json");
    diff_objectives_json(a / "objectives.json", b / "objectives.json");
    if(exists(a / "normalization.json") || exists(b / "normalization.json"))
      {
        if(check_normalization)
          diff_normalization_json(a / "normalization.json",
                                  b / "normalization.json");
      }
    Parse_Control_Json control(a / "control.json");
    for(int i = 0; i < control.num_blocks; ++i)
      {
        auto block_info_name = "block_info_" + std::to_string(i) + ".json";
        diff_block_info_json(a / block_info_name, b / block_info_name);

        auto block_data_name = "block_data_" + std::to_string(i);
        diff_block_data_file(a / block_data_name, b / block_data_name);
      }
  }
}
