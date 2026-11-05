#pragma once

#include "Abstract_Json_Element_Parser.hxx"

#include <functional>

class Json_Int64_Parser final : public Abstract_Json_Element_Parser<int64_t>
{
  int64_t result = 0;

public:
  using Abstract_Json_Element_Parser::value_type;

  Json_Int64_Parser(
    const bool skip, const std::function<void(value_type &&)> &on_parsed,
    const std::function<void()> &on_skipped = [] {})
      : Abstract_Json_Element_Parser(skip, on_parsed, on_skipped)
  {}

  bool json_int64(const int64_t value) override
  {
    if(!this->skip)
      this->result = value;

    this->on_end();
    return true;
  }

  bool json_uint64(const uint64_t value) override { return json_int64(value); }
  bool json_int(const int value) override { return json_int64(value); }
  bool json_uint(const unsigned value) override { return json_int64(value); }

  void reset(const bool skip) override { this->skip = skip; }
  value_type get_result() override { return result; }
};
