#pragma once

#include "sdpb_util/Number_State.hxx"
#include "sdpb_util/Vector_State.hxx"

#include <El.hpp>
#include <rapidjson/reader.h>

using namespace std::string_literals;
struct Checkpoint_Parser
    : public rapidjson::BaseReaderHandler<rapidjson::UTF8<>, Checkpoint_Parser>
{
  bool inside = false, parsing_generation = false, parsing_threshold = false,
       parsing_c_scale = false, parsing_yp = false, parsing_y_transform = false,
       parsing_b = false, parsing_points = false;

  Number_State<El::BigFloat> generation_state, threshold_state, c_scale_state;
  Vector_State<Number_State<El::BigFloat>> yp_state, b_state;
  Vector_State<Vector_State<Number_State<El::BigFloat>>> points_state,
    y_transform_state;

  Checkpoint_Parser()
      : generation_state("generation"s), threshold_state("threshold"s),
        c_scale_state("c_scale"s), yp_state({"yp"s, ""s}), b_state({"b"s, ""s}),
        points_state({"points"s, ""s, ""s}),
        y_transform_state({"y_transform"s, ""s, ""s})
  {}

  bool Null() { throw std::runtime_error("Null not allowed"); }
  bool Bool(bool) { throw std::runtime_error("Bool not allowed"); }

  bool Int(int)
  {
    throw std::runtime_error(
      "Int not allowed.  You must quote all numbers as strings.");
  }
  bool Uint(unsigned)
  {
    throw std::runtime_error(
      "Uint not allowed.  You must quote all numbers as strings.");
  }
  bool Int64(int64_t)
  {
    throw std::runtime_error(
      "Int64 not allowed.  You must quote all numbers as strings.");
  }
  bool Uint64(uint64_t)
  {
    throw std::runtime_error(
      "Uint64 not allowed.  You must quote all numbers as strings.");
  }
  bool Double(double)
  {
    throw std::runtime_error(
      "Double not allowed.  You must quote all numbers as strings.");
  }
  bool RawNumber(const Ch *, rapidjson::SizeType, bool)
  {
    throw std::runtime_error(
      "Numbers not allowed.  You must quote all numbers as strings.");
  }
  bool String(const Ch *str, rapidjson::SizeType length, bool copy);
  bool StartObject();
  bool Key(const Ch *str, rapidjson::SizeType length, bool copy);
  bool EndObject(rapidjson::SizeType memberCount);
  bool StartArray();
  bool EndArray(rapidjson::SizeType elementCount);
};
