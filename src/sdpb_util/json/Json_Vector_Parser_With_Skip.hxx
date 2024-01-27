#pragma once

#include "Abstract_Json_Array_Parser_With_Skip.hxx"
#include "Vector_Parse_Result_With_Skip.hxx"

// Parses JSON array elements with TElementParser,
// skipping them according to skip_element(index) == true,
// and collects them to Vector_Parse_Result_With_Skip<TElementParser::value_type>
template <class TElementParser>
class Json_Vector_Parser_With_Skip final
    : public Abstract_Json_Array_Parser_With_Skip<
        Vector_Parse_Result_With_Skip<typename TElementParser::value_type>,
        TElementParser>
{
public:
  ~Json_Vector_Parser_With_Skip() = default;

public:
  using element_type = typename TElementParser::value_type;
  using value_type = Vector_Parse_Result_With_Skip<element_type>;

private:
  value_type result;

public:
  template <class... TArgs>
  Json_Vector_Parser_With_Skip(
    const bool skip,
    const std::function<void(Vector_Parse_Result_With_Skip<element_type> &&)>
      &on_parsed,
    const std::function<void()> &on_skipped,
    const std::function<bool(size_t index)> &skip_element,
    const TArgs &...element_parser_args)
      : Abstract_Json_Array_Parser_With_Skip<
          Vector_Parse_Result_With_Skip<element_type>, TElementParser>(
          skip, on_parsed, on_skipped, skip_element, element_parser_args...)
  {}
  void on_element_parsed(element_type &&value, size_t index) override
  {
    this->result.add_element(std::move(value), index);
  }
  void on_element_skipped(size_t index) override
  {
    this->result.skip_element(index);
  }
  value_type get_result() override { return std::move(result); }
  void clear_result() override { result.reset(); }
};
