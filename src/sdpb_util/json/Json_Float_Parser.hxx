#pragma once

#include "Abstract_Element_Parser.hxx"
#include "Json_Vector_Parser.hxx"

#include <El.hpp>

#include <rapidjson/reader.h>
#include <rapidjson/encodings.h>

#include <functional>

// TODO merge with NumberState
template <class TFloat>
class Json_Float_Parser final : public Abstract_Json_Element_Parser<TFloat>
{
  using SizeType = rapidjson::SizeType;
  using Ch = rapidjson::UTF8<>::Ch;

  TFloat result;

public:
  using typename Abstract_Json_Element_Parser<TFloat>::value_type;

  Json_Float_Parser(
    const bool skip, const std::function<void(value_type &&)> &on_parsed,
    const std::function<void()> &on_skipped = [] {})
      : Abstract_Json_Element_Parser<TFloat>(skip, on_parsed, on_skipped)
  {}

  bool json_string(const rapidjson::UTF8<>::Ch *str,
                   rapidjson::SizeType length, bool /*copy*/) override
  {
    if(!this->skip)
      {
        std::string string_value(str, length);
        try
          {
            // GMP does not have inf or nan, so we approximate inf
            // with max double.
            // TODO throw an error insted?
            using namespace std::string_literals;
            this->result = string_value == "inf"s
                             ? TFloat(std::numeric_limits<double>::max())
                             : TFloat(string_value);
          }
        catch(std::exception &e)
          {
            RUNTIME_ERROR("Json_Number_Parser failed to parse '", string_value,
                          "': ", e.what());
          }
        catch(...)
          {
            RUNTIME_ERROR("Json_Number_Parser failed to parse '", string_value,
                          "'");
          }
      }

    this->on_end();
    return true;
  }

  void reset(const bool skip) override { this->skip = skip; }
  value_type get_result() override { return std::move(result); }
};

template <class TFloat>
using Json_Float_Vector_Parser = Json_Vector_Parser<Json_Float_Parser<TFloat>>;
