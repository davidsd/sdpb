#pragma once

#include "Float.hxx"
#include "Test_Case_Runner.hxx"

#include <catch2/catch_amalgamated.hpp>
#include <boost/filesystem.hpp>
#include <rapidjson/document.h>

namespace Test_Util::Json
{
  using Json_Value = rapidjson::GenericValue<rapidjson::UTF8<>>;

  inline Float parse_Float(const Json_Value &json_value)
  {
    return Float(json_value.GetString());
  }
  inline Float_Vector parse_Float_Vector(const Json_Value &json_value)
  {
    Float_Vector result;
    for(const auto &element : json_value.GetArray())
      {
        result.emplace_back(parse_Float(element));
      }
    return result;
  }
  inline Float_Matrix parse_Float_Matrix(const Json_Value &json_value)
  {
    Float_Matrix result;
    for(const auto &element : json_value.GetArray())
      {
        result.emplace_back(parse_Float_Vector(element));
      }
    return result;
  }
}
