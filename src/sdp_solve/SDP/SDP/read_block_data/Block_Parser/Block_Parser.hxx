#pragma once

#include "sdpb_util/Vector_State.hxx"
#include "sdpb_util/Number_State.hxx"

#include <rapidjson/reader.h>

#include <string>

using namespace std::string_literals;
struct Block_Parser
    : public rapidjson::BaseReaderHandler<rapidjson::UTF8<>, Block_Parser>
{
  bool inside = false, parsing_c = false, parsing_B = false,
       parsing_bilinear_bases_even = false, parsing_bilinear_bases_odd = false;

  Vector_State<Number_State<El::BigFloat>> c_state;
  Vector_State<Vector_State<Number_State<El::BigFloat>>> B_state,
    bilinear_bases_even_state, bilinear_bases_odd_state;

  Block_Parser()
      : c_state({"c"s, ""s}),
        B_state({"B"s, ""s, ""s}),
        bilinear_bases_even_state({"bilinear_bases_even"s, ""s, ""s}),
        bilinear_bases_odd_state({"bilinear_bases_odd"s, ""s, ""s})
  {}

  bool Null() { RUNTIME_ERROR("Null not allowed"); }
  bool Bool(bool) { RUNTIME_ERROR("Bool not allowed"); }
  template <typename T> bool parse_integer(const T &i)
  {
    RUNTIME_ERROR("Integer not allowed");
  }
  bool Int(int i) { return parse_integer(i); }
  bool Uint(unsigned i) { return parse_integer(i); }
  bool Int64(int64_t i) { return parse_integer(i); }
  bool Uint64(uint64_t i) { return parse_integer(i); }

  bool Double(double d)
  {
    RUNTIME_ERROR("Invalid input '", d,
                  "'. All floating point numbers must be quoted as strings.");
  }
  bool RawNumber(const Ch *characters, rapidjson::SizeType size, bool)
  {
    RUNTIME_ERROR("Invalid input '", std::string(characters, size),
                  "'.  All floating point numbers must be quoted as strings.");
  }
  bool StartObject()
  {
    if(!inside)
      {
        inside = true;
      }
    else
      {
        RUNTIME_ERROR(
          "Invalid input. No objects allowed inside the main object.");
      }
    return true;
  }
  bool EndObject(rapidjson::SizeType)
  {
    inside = false;
    return true;
  }

  bool String(const Ch *str, rapidjson::SizeType length, bool copy);
  bool Key(const Ch *str, rapidjson::SizeType length, bool copy);
  bool StartArray();
  bool EndArray(rapidjson::SizeType elementCount);
};
