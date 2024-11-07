#include "Parse_SDP.hxx"

#include "diff.hxx"
#include "json.hxx"
#include "sdpb_util/boost_serialization.hxx"

#include <fstream>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>

namespace fs = std::filesystem;

using namespace Test_Util::Json;

Parse_Control_Json::Parse_Control_Json(const fs::path &path)
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

Parse_Pmp_Info_Json::Parse_Pmp_Info_Json(const fs::path &path)
{
  CAPTURE(path);
  REQUIRE(exists(path));
  std::ifstream is(path);
  rapidjson::IStreamWrapper wrapper(is);
  rapidjson::Document document;
  document.ParseStream(wrapper);

  REQUIRE(document.IsArray());
  for(const auto &item : document.GetArray())
    {
      auto &pvm_info = pmp_info.emplace_back();
      pvm_info.block_path = item["block_path"].GetString();
      pvm_info.prefactor = parse_Damped_Rational(item["prefactor"]);
      pvm_info.reduced_prefactor
        = parse_Damped_Rational(item["reducedPrefactor"]);
      pvm_info.sample_points = parse_Float_Vector(item["samplePoints"]);
      pvm_info.sample_scalings = parse_Float_Vector(item["sampleScalings"]);
      pvm_info.reduced_sample_scalings
        = parse_Float_Vector(item["reducedSampleScalings"]);
    }
}

Parse_Objectives_Json::Parse_Objectives_Json(const fs::path &path)
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

Parse_Normalization_Json::Parse_Normalization_Json(const fs::path &path)
{
  CAPTURE(path);
  REQUIRE(exists(path));
  std::ifstream is(path);
  rapidjson::IStreamWrapper wrapper(is);
  rapidjson::Document document;
  document.ParseStream(wrapper);

  normalization = parse_Float_Vector(document["normalization"]);
}

Parse_Block_Info_Json::Parse_Block_Info_Json(const fs::path &path)
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

Parse_Block_Data::Parse_Block_Data(const fs::path &block_path_no_extension)
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

void Parse_Block_Data::parse_json(const fs::path &block_path)
{
  CAPTURE(block_path);
  std::ifstream is(block_path);
  rapidjson::IStreamWrapper wrapper(is);
  rapidjson::Document document;
  document.ParseStream(wrapper);

  bilinear_bases_even = parse_Float_Matrix(document["bilinear_bases_even"]);
  bilinear_bases_odd = parse_Float_Matrix(document["bilinear_bases_odd"]);
  constraint_constants = parse_Float_Vector(document["c"]);
  constraint_matrix = parse_Float_Matrix(document["B"]);
}
void Parse_Block_Data::parse_bin(const fs::path &block_path)
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

Parse_SDP::Parse_SDP(const fs::path &sdp_path,
                     const Test_Util::Test_Case_Runner &runner)
    : sdp_dir(is_regular_file(sdp_path) ? runner.unzip_to_temp_dir(sdp_path)
                                        : sdp_path),
      control(sdp_dir / "control.json"),
      pmp_info(sdp_dir / "pmp_info.json"),
      objectives(sdp_dir / "objectives.json")
{
  if(exists(sdp_dir / "normalization.json"))
    normalization.emplace(sdp_dir / "normalization.json");
  for(size_t i = 0; i < control.num_blocks; ++i)
    {
      block_info.emplace_back(sdp_dir
                              / ("block_info_" + std::to_string(i) + ".json"));
      block_data.emplace_back(sdp_dir / ("block_data_" + std::to_string(i)));
    }
}