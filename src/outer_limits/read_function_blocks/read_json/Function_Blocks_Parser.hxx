#pragma once

#include "Function_State.hxx"

#include <El.hpp>
#include <rapidjson/reader.h>

using namespace std::string_literals;
struct Function_Blocks_Parser
    : public rapidjson::BaseReaderHandler<rapidjson::UTF8<>,
                                          Function_Blocks_Parser>
{
  bool inside = false, parsing_objective = false,
       parsing_normalization = false, parsing_functions = false;

  Vector_State<Number_State<El::BigFloat>> objective_state,
    normalization_state;
  Vector_State<Vector_State<Vector_State<Vector_State<Function_State>>>>
    functions_state;

  Function_Blocks_Parser()
      : objective_state({"objective"s, ""s}),
        normalization_state({"normalization"s, ""s}),
        functions_state({"functions"s, ""s, ""s, ""s, ""s})
  {}

  bool Null() { RUNTIME_ERROR("Null not allowed"); }
  bool Bool(bool) { RUNTIME_ERROR("Bool not allowed"); }
  bool Int(int)
  {
    RUNTIME_ERROR(
      "Int not allowed.  You must quote all numbers as strings.");
  }
  bool Uint(unsigned)
  {
    RUNTIME_ERROR(
      "Uint not allowed.  You must quote all numbers as strings.");
  }
  bool Int64(int64_t)
  {
    RUNTIME_ERROR(
      "Int64 not allowed.  You must quote all numbers as strings.");
  }
  bool Uint64(uint64_t)
  {
    RUNTIME_ERROR(
      "Uint64 not allowed.  You must quote all numbers as strings.");
  }
  bool Double(double)
  {
    RUNTIME_ERROR(
      "Double not allowed.  You must quote all numbers as strings.");
  }
  bool RawNumber(const Ch *, rapidjson::SizeType, bool)
  {
    RUNTIME_ERROR(
      "Numbers not allowed.  You must quote all numbers as strings.");
  }
  bool String(const Ch *str, rapidjson::SizeType length, bool copy);
  bool StartObject();
  bool Key(const Ch *str, rapidjson::SizeType length, bool copy);
  bool EndObject(rapidjson::SizeType memberCount);
  bool StartArray();
  bool EndArray(rapidjson::SizeType elementCount);
};
