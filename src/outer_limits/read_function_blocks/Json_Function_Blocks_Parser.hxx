#pragma once

#include "Functions_File_Parse_Result.hxx"
#include "Json_Function_Parser.hxx"
#include "outer_limits/Function.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/json/Abstract_Json_Object_Parser.hxx"
#include "sdpb_util/json/Json_Float_Parser.hxx"
#include "sdpb_util/Boost_Float.hxx"

using namespace std::string_literals;
struct Json_Functions_Parser final
    : Abstract_Json_Object_Parser<Functions_File_Parse_Result>
{
private:
  Functions_File_Parse_Result result;
  Json_Float_Vector_Parser<El::BigFloat> objective_parser;
  Json_Float_Vector_Parser<El::BigFloat> normalization_parser;
  Json_Vector_Parser<Json_Vector_Parser<
    Json_Vector_Parser<Json_Vector_Parser<Json_Function_Parser>>>>
    functions_parser;

  // Abstract_Json_Object_Parser implementation
public:
  Functions_File_Parse_Result get_result() override
  {
    return std::move(result);
  }

protected:
  Abstract_Json_Reader_Handler &element_parser(const std::string &key) override
  {
    if(key == "objective")
      return objective_parser;
    if(key == "normalization")
      return normalization_parser;
    if(key == "functions")
      return functions_parser;
    RUNTIME_ERROR("Unknown key: '", key, "'");
  }

public:
  void reset_element_parsers(bool skip) override
  {
    objective_parser.reset(skip);
    normalization_parser.reset(skip);
    functions_parser.reset(skip);
  }
  void clear_result() override
  {
    result.objective.clear();
    result.normalization.clear();
    result.functions.clear();
  }

public:
  [[nodiscard]] Json_Functions_Parser(
    const std::function<void(Functions_File_Parse_Result &&)> &on_parsed)
      : Abstract_Json_Object_Parser(false, on_parsed, [] {}),
        objective_parser(false,
                         [this](std::vector<El::BigFloat> &&objective) {
                           this->result.objective = std::move(objective);
                         }),
        normalization_parser(
          false,
          [this](std::vector<El::BigFloat> &&normalization) {
            this->result.normalization = std::move(normalization);
          }),
        functions_parser(
          false,
          [this](std::vector<std::vector<std::vector<std::vector<Function>>>>
                   &&functions) {
            this->result.functions = std::move(functions);
          })
  {}
};
