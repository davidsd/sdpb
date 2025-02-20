#pragma once

#include "Abstract_Json_Element_Parser.hxx"
#include "sdpb_util/assert.hxx"

#include <rapidjson/encodings.h>

#include <functional>

template <class TResult>
class Json_String_Element_Parser : public Abstract_Json_Element_Parser<TResult>
{
  using SizeType = rapidjson::SizeType;
  using Ch = rapidjson::UTF8<>::Ch;

  TResult result;

protected:
  virtual TResult from_string(const std::string &string_value)
  {
    return TResult(string_value);
  }

public:
  using typename Abstract_Json_Element_Parser<TResult>::value_type;

  virtual ~Json_String_Element_Parser() = default;

  Json_String_Element_Parser(
    const bool skip, const std::function<void(value_type &&)> &on_parsed,
    const std::function<void()> &on_skipped = [] {})
      : Abstract_Json_Element_Parser<TResult>(skip, on_parsed, on_skipped)
  {}

  bool json_string(const rapidjson::UTF8<>::Ch *str,
                   rapidjson::SizeType length, bool /*copy*/) override
  {
    if(!this->skip)
      {
        std::string string_value(str, length);
        try
          {
            this->result = from_string(string_value);
          }
        catch(std::exception &e)
          {
            RUNTIME_ERROR("'", typeid(*this).name(),
                          "' failed to parse string \"", string_value,
                          "\": ", e.what());
          }
        catch(...)
          {
            RUNTIME_ERROR("'", typeid(*this).name(),
                          "' failed to parse string \"", string_value, "\"");
          }
      }

    this->on_end();
    return true;
  }

  void reset(const bool skip) override { this->skip = skip; }
  value_type get_result() override { return std::move(result); }
};
