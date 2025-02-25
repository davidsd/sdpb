#pragma once

#define RAPIDJSON_ASSERT(x) ASSERT(x)
#define RAPIDJSON_NOEXCEPT_ASSERT(x) assert(x)

#include "Float.hxx"
#include "pmp/Damped_Rational.hxx"
#include "sdpb_util/assert.hxx"

#include <catch2/catch_amalgamated.hpp>
#include <filesystem>
#include <rapidjson/document.h>

namespace Test_Util::Json
{
  using Json_Value = rapidjson::GenericValue<rapidjson::UTF8<>>;

  template <class TValue> TValue parse(const Json_Value &json_value)
  {
    return TValue(json_value.GetString());
  }

  template <class TValue>
  std::vector<TValue>
  parse_vector(const Json_Value &json_value,
               const std::function<TValue(const Json_Value &)> parse_element
               = parse<TValue>)
  {
    std::vector<TValue> result;
    for(const auto &element : json_value.GetArray())
      {
        result.emplace_back(parse_element(element));
      }
    return result;
  }

  inline Float parse_Float(const Json_Value &json_value)
  {
    return parse<Float>(json_value);
  }
  inline Float_Vector parse_Float_Vector(const Json_Value &json_value)
  {
    return parse_vector<Float>(json_value);
  }
  inline Float_Matrix parse_Float_Matrix(const Json_Value &json_value)
  {
    auto rows_array = json_value.GetArray();
    int height = rows_array.Size();
    int width = height == 0 ? 1 : rows_array[0].GetArray().Size();
    Float_Matrix result(height, width);

    for(int row = 0; row < height; ++row)
      {
        auto cols_array = rows_array[row].GetArray();
        for(int col = 0; col < width; ++col)
          {
            Float value = parse_Float(cols_array[col]);
            result.Set(row, col, value);
          }
      }
    return result;
  }

  inline Damped_Rational parse_Damped_Rational(const Json_Value &json_value)
  {
    Damped_Rational result;
    auto obj = json_value.GetObject();
    result.base = parse<Boost_Float>(obj["base"]);
    result.constant = parse<Boost_Float>(obj["constant"]);
    result.poles = parse_vector<Boost_Float>(obj["poles"]);
    return result;
  }
}
