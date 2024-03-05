#pragma once

#include "Abstract_Json_Array_Parser_With_Skip.hxx"

#include <vector>

// Parses JSON array elements with TElementParser
// and collects them to TResult
template <class TResult, class TElementParser>
class Abstract_Json_Vector_Parser
    : public Abstract_Json_Array_Parser_With_Skip<TResult, TElementParser>
{
protected:
  ~Abstract_Json_Vector_Parser() = default;

public:
  using element_type = typename TElementParser::value_type;
  using value_type = TResult;

  template <class... TArgs>
  Abstract_Json_Vector_Parser(
    const bool skip, const std::function<void(value_type &&)> &on_parsed,
    const std::function<void()> &on_skipped,
    const TArgs &...element_parser_args)
      : Abstract_Json_Array_Parser_With_Skip<value_type, TElementParser>(
          skip, on_parsed, on_skipped,
          // Never skip individual elements:
          [](size_t /*index*/) { return false; },
          element_parser_args...)
  {}
  void on_element_skipped(size_t /*index*/) override {}
};

// Parses JSON array elements with TElementParser
// and collects them to std::vector<TElementParser::value_type>
template <class TElementParser>
class Json_Vector_Parser final
    : public Abstract_Json_Vector_Parser<
        std::vector<typename TElementParser::value_type>, TElementParser>
{
public:
  virtual ~Json_Vector_Parser() = default;
  using base_type = Abstract_Json_Vector_Parser<
    std::vector<typename TElementParser::value_type>, TElementParser>;
  using element_type = typename base_type::element_type;
  using value_type = typename base_type::value_type;

private:
  value_type result;

public:
  Json_Vector_Parser(
    const bool skip, const std::function<void(value_type &&)> &on_parsed,
    const std::function<void()> &on_skipped = [] {})
      : base_type(skip, on_parsed, on_skipped)
  {}

  void on_element_parsed(element_type &&value, size_t index) override
  {
    ASSERT_EQUAL(index , this->result.size());
    this->result.push_back(std::move(value));
  }

private:
  void clear_result() override { result.clear(); }
  value_type get_result() override { return std::move(result); }
};
