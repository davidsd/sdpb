#pragma once

#include "Positive_Matrix_With_Prefactor_State.hxx"

#include <rapidjson/reader.h>

using namespace std::string_literals;
struct JSON_Parser
    : public rapidjson::BaseReaderHandler<rapidjson::UTF8<>, JSON_Parser>
{
  bool inside = false, parsing_objective = false,
       parsing_normalization = false,
       parsing_positive_matrices_with_prefactor = false;

  Vector_State<Number_State<El::BigFloat>> objective_state,
    normalization_state;
  Vector_State<Positive_Matrix_With_Prefactor_State>
    positive_matrices_with_prefactor_state;

  JSON_Parser()
      : objective_state({"objective"s, ""s}),
        normalization_state({"normalization"s, ""s}),
        positive_matrices_with_prefactor_state(
          {"PositiveMatrixWithPrefactorArray"s, ""s, "DampedRational"s,
           "polynomials"s, ""s, ""s, ""s, ""s})
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
