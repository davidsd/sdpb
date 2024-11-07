#pragma once
#include "Float.hxx"
#include "Test_Case_Runner.hxx"
#include "pmp/PVM_Info.hxx"

#include <optional>
#include <vector>
#include <boost/core/noncopyable.hpp>
#include <catch2/catch_amalgamated.hpp>
#include <filesystem>

// Parser classes
// TODO replace them with faster parsers from SDPB?
// These parsers use DOM, which is simpler but slower in comparison to SAX

struct Parse_Control_Json : boost::noncopyable
{
  int num_blocks;
  std::string command;
  explicit Parse_Control_Json(const std::filesystem::path &path);
};

struct Parse_Pmp_Info_Json : boost::noncopyable
{
  std::vector<PVM_Info> pmp_info;
  explicit Parse_Pmp_Info_Json(const std::filesystem::path &path);
};

struct Parse_Objectives_Json : boost::noncopyable
{
  Float constant;
  std::vector<Float> b;
  explicit Parse_Objectives_Json(const std::filesystem::path &path);
};

struct Parse_Normalization_Json : boost::noncopyable
{
  std::vector<Float> normalization;
  explicit Parse_Normalization_Json(const std::filesystem::path &path);
};

struct Parse_Block_Info_Json
{
  size_t dim;
  size_t num_points;
  explicit Parse_Block_Info_Json(const std::filesystem::path &path);
};

struct Parse_Block_Data
{
  Float_Matrix bilinear_bases_even;
  Float_Matrix bilinear_bases_odd;
  Float_Vector constraint_constants;
  Float_Matrix constraint_matrix;
  explicit
  Parse_Block_Data(const std::filesystem::path &block_path_no_extension);

private:
  void parse_json(const std::filesystem::path &block_path);
  void parse_bin(const std::filesystem::path &block_path);
};

struct Parse_SDP
{
  const std::filesystem::path sdp_dir;
  Parse_Control_Json control;
  Parse_Pmp_Info_Json pmp_info;
  Parse_Objectives_Json objectives;
  std::optional<Parse_Normalization_Json> normalization;
  std::vector<Parse_Block_Info_Json> block_info;
  std::vector<Parse_Block_Data> block_data;

  Parse_SDP(const std::filesystem::path &sdp_path,
            const Test_Util::Test_Case_Runner &runner);
};