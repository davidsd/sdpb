#pragma once

#include "Json_Positive_Matrix_With_Prefactor_Parser.hxx"
#include "pmp_read/PMP_File_Parse_Result.hxx"
#include "sdpb_util/json/Abstract_Json_Object_Parser.hxx"
#include "sdpb_util/json/Json_Float_Parser.hxx"
#include "sdpb_util/json/Json_Vector_Parser_With_Skip.hxx"

class Json_PMP_Parser final
    : public Abstract_Json_Object_Parser<PMP_File_Parse_Result>
{
  using SizeType = rapidjson::SizeType;
  using Ch = rapidjson::UTF8<>::Ch;

private:
  using BigFloat_Vector_Parser = Json_Vector_Parser<Json_BigFloat_Parser>;

  using Json_Positive_Matrix_With_Prefactor_Array_Parser
    = Json_Vector_Parser_With_Skip<Json_Positive_Matrix_With_Prefactor_Parser>;

  PMP_File_Parse_Result result;

  // Nested parsers
  BigFloat_Vector_Parser objective_parser;
  BigFloat_Vector_Parser normalization_parser;
  Json_Positive_Matrix_With_Prefactor_Array_Parser matrices_parser;

public:
  Json_PMP_Parser(
    bool should_parse_objective, bool should_parse_normalization,
    const std::function<bool(size_t matrix_index)> &should_parse_matrix,
    const std::function<void(PMP_File_Parse_Result &&result)> &on_parsed);
  Abstract_Json_Reader_Handler &
  element_parser(const std::string &key) override;

public:
  void clear_result() override;
  value_type get_result() override;
  void reset_element_parsers(bool skip) override;
};
