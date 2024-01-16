#pragma once

#include "Abstract_Json_Element_Parser.hxx"

#include <functional>

// TODO merge with VectorState
template <class TResult, class TElementParser>
class Abstract_Json_Array_Parser_With_Skip
    : public Abstract_Json_Element_Parser<TResult>
{
public:
  using element_type = typename TElementParser::value_type;
  using value_type = TResult;
  using base_type = Abstract_Json_Element_Parser<TResult>;

protected:
  ~Abstract_Json_Array_Parser_With_Skip() = default;

private:
  using SizeType = rapidjson::SizeType;
  using Ch = rapidjson::UTF8<>::Ch;

  std::function<bool(size_t index)> skip_element;

  enum State
  {
    Start,
    Inside,
    // TODO maybe remove End? We already have on_end() callbacks.
    End
  } state
    = Start;

  size_t index = 0;
  int array_level = 0;

  TElementParser element_parser;

public:
  template <class... TArgs>
  Abstract_Json_Array_Parser_With_Skip(
    bool skip, const std::function<void(TResult &&)> &on_parsed,
    const std::function<void()> &on_skipped,
    const std::function<bool(size_t index)> &skip_element,
    const TArgs &...element_parser_args);

private:
  virtual void on_element_parsed(element_type &&value, size_t index) = 0;
  virtual void on_element_skipped(size_t index) = 0;

public:
  void reset(const bool skip) override
  {
    this->skip = skip;
    element_parser.reset(this->skip_element(0));
    state = Start;
    index = 0;
    array_level = 0;
    clear_result();
  }

protected:
  virtual void clear_result() = 0;

#define ABSTRACT_JSON_ARRAY_ELEMENT_PARSER(func)                              \
  (state == Inside ? element_parser.func : base_type::func)

public:
  bool json_null() override
  {
    return ABSTRACT_JSON_ARRAY_ELEMENT_PARSER(json_null());
  }
  bool json_bool(bool b) override
  {
    return ABSTRACT_JSON_ARRAY_ELEMENT_PARSER(json_bool(b));
  }
  bool json_int(int i) override
  {
    return ABSTRACT_JSON_ARRAY_ELEMENT_PARSER(json_int(i));
  }
  bool json_uint(unsigned i) override
  {
    return ABSTRACT_JSON_ARRAY_ELEMENT_PARSER(json_uint(i));
  }
  bool json_int64(int64_t i) override
  {
    return ABSTRACT_JSON_ARRAY_ELEMENT_PARSER(json_int64(i));
  }
  bool json_uint64(uint64_t i) override
  {
    return ABSTRACT_JSON_ARRAY_ELEMENT_PARSER(json_uint64(i));
  }
  bool json_double(double d) override
  {
    return ABSTRACT_JSON_ARRAY_ELEMENT_PARSER(json_double(d));
  }
  bool json_raw_number(const Ch *str, SizeType length, bool copy) override
  {
    return ABSTRACT_JSON_ARRAY_ELEMENT_PARSER(
      json_raw_number(str, length, copy));
  }
  bool json_string(const Ch *str, SizeType length, bool copy) override
  {
    return ABSTRACT_JSON_ARRAY_ELEMENT_PARSER(json_string(str, length, copy));
  }
  bool json_start_object() override
  {
    return ABSTRACT_JSON_ARRAY_ELEMENT_PARSER(json_start_object());
  }
  bool json_key(const Ch *str, SizeType length, bool copy) override
  {
    return ABSTRACT_JSON_ARRAY_ELEMENT_PARSER(json_key(str, length, copy));
  }
  bool json_end_object(SizeType memberCount) override
  {
    return ABSTRACT_JSON_ARRAY_ELEMENT_PARSER(json_end_object(memberCount));
  }
  bool json_start_array() override
  {
    ++array_level;
    switch(state)
      {
      case Start:
        // NB: we set element_parser.skip in constructor,
        // so this is not necessary here.
        // TODO: remove skip from constructors (set false by default)?
        // This would make code less verbose,
        // but feels somewhat less safe.
        element_parser.reset(this->skip_element(0));
        state = Inside;
        return true;
      case Inside:
        return ABSTRACT_JSON_ARRAY_ELEMENT_PARSER(json_start_array());
      default: return base_type::json_start_array();
      }
  }
  bool json_end_array(SizeType elementCount) override
  {
    --array_level;
    switch(state)
      {
        case Inside: {
          if(array_level == 0)
            {
              state = End;
              // TODO callbacks in on_end() can change state!
              // e.g. it always happens in case of nested arrays
              // Feels like spaghetti, need refactoring
              this->on_end();
              return true;
            }
          return ABSTRACT_JSON_ARRAY_ELEMENT_PARSER(
            json_end_array(elementCount));
        }
      default: return base_type::json_start_array();
      }
  }
};
template <class TResult, class TElementParser>
template <class... TArgs>
Abstract_Json_Array_Parser_With_Skip<TResult, TElementParser>::
  Abstract_Json_Array_Parser_With_Skip(
    bool skip, const std::function<void(TResult &&)> &on_parsed,
    const std::function<void()> &on_skipped,
    const std::function<bool(size_t index)> &skip_element_func,
    const TArgs &...element_parser_args)
    : base_type(skip, on_parsed, on_skipped),
      skip_element([this, skip_element_func](size_t index) {
        return this->skip || skip_element_func(index);
      }),
      element_parser(
        this->skip_element(0),
        [this](element_type &&value) {
          this->on_element_parsed(std::forward<element_type>(value),
                                  this->index);
          ++this->index;
          this->element_parser.reset(this->skip_element(this->index));
        },
        [this] {
          this->on_element_skipped(this->index);
          ++this->index;
          this->element_parser.reset(this->skip_element(this->index));
        },
        element_parser_args...)
{}
