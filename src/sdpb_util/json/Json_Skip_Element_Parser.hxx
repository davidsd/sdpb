#pragma once

#include "Abstract_Json_Element_Parser.hxx"

// Dummy parser that simply skips everything.
// Conceptually, it should be Abstract_Json_Element_Parser<void>,
// but we cannot use void as a template argument,
// so we use bool and return true (which does not mean anything).
class Json_Skip_Element_Parser final
    : public Abstract_Json_Element_Parser<bool>
{
public:
  using typename Abstract_Json_Element_Parser::value_type;

  explicit Json_Skip_Element_Parser(
    const bool skip = true,
    const std::function<void(value_type &&)> &on_parsed
    = [](value_type &&value) { return value; },
    const std::function<void()> &on_skipped = [] {})
      : Abstract_Json_Element_Parser(true, on_parsed, on_skipped)
  {}

  bool json_default() override { return true; }
  bool json_null() override { return true; }
  bool json_bool(bool) override { return true; }
  bool json_int(int) override { return true; }
  bool json_uint(unsigned) override { return true; }
  bool json_int64(int64_t) override { return true; }
  bool json_uint64(uint64_t) override { return true; }
  bool json_double(double) override { return true; }
  bool json_raw_number(const Ch *, SizeType, bool) override { return true; }
  bool json_string(const Ch *, SizeType, bool) override { return true; }
  bool json_start_object() override { return true; }
  bool json_key(const Ch *, SizeType, bool) override { return true; }
  bool json_end_object(SizeType) override { return true; }
  bool json_start_array() override { return true; }
  bool json_end_array(SizeType) override { return true; }

  void reset(const bool skip) override { this->skip = skip; }
  value_type get_result() override { return true; }
};
