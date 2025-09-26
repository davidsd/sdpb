#pragma once

#include "Abstract_Json_Element_Parser.hxx"

#include <El.hpp>

#include <rapidjson/reader.h>
#include <rapidjson/encodings.h>

#include <functional>
#include <string>

class Json_String_Parser final : public Abstract_Json_Element_Parser<std::string>
{
  using SizeType = rapidjson::SizeType;
  using Ch = rapidjson::UTF8<>::Ch;

  std::string result;

public:
  using typename Abstract_Json_Element_Parser::value_type;

  Json_String_Parser(
    const bool skip, const std::function<void(value_type &&)> &on_parsed,
    const std::function<void()> &on_skipped = [] {})
      : Abstract_Json_Element_Parser(skip, on_parsed, on_skipped)
  {}

  bool json_string(const rapidjson::UTF8<>::Ch *str,
                   rapidjson::SizeType length, bool /*copy*/) override
  {
    if(!this->skip)
      this->result = std::string(str, length);

    this->on_end();
    return true;
  }

  void reset(const bool skip) override { this->skip = skip; }
  value_type get_result() override { return std::move(result); }
};

