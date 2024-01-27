#pragma once

#include "Abstract_Json_Element_Parser.hxx"
#include "Abstract_Json_Reader_Handler.hxx"

#include <functional>

template <class TResult>
class Abstract_Json_Object_Parser
    : public Abstract_Json_Element_Parser<TResult>
{
public:
  using value_type = TResult;

protected:
  ~Abstract_Json_Object_Parser() = default;

private:
  using SizeType = rapidjson::SizeType;
  using Encoding = rapidjson::UTF8<>;
  using Ch = rapidjson::UTF8<>::Ch;
  using base_type = Abstract_Json_Element_Parser<TResult>;

  enum State
  {
    Start,
    Inside,
    // InsideElement,
    End
  } state
    = Start;

  int object_level = 0;
  int array_level = 0;
  std::string key;

protected:
  virtual Abstract_Json_Reader_Handler &element_parser(const std::string &key)
    = 0;

public:
  Abstract_Json_Object_Parser(
    bool skip, const std::function<void(TResult &&)> &on_parsed,
    const std::function<void()> &on_skipped = [] {});

  void reset(bool skip) override;
  virtual void reset_element_parsers(bool skip) = 0;
  virtual void clear_result() = 0;

#define ABSTRACT_JSON_OBJECT_ELEMENT_PARSER(member)                           \
  (state == Inside ? element_parser(key).member : base_type::member)

public:
  bool json_default() override
  {
    return ABSTRACT_JSON_OBJECT_ELEMENT_PARSER(json_default());
  }
  bool json_null() override
  {
    return ABSTRACT_JSON_OBJECT_ELEMENT_PARSER(json_null());
  }
  bool json_bool(bool b) override
  {
    return ABSTRACT_JSON_OBJECT_ELEMENT_PARSER(json_bool(b));
  }
  bool json_int(int i) override
  {
    return ABSTRACT_JSON_OBJECT_ELEMENT_PARSER(json_int(i));
  }
  bool json_uint(unsigned i) override
  {
    return ABSTRACT_JSON_OBJECT_ELEMENT_PARSER(json_uint(i));
  }
  bool json_int64(int64_t i) override
  {
    return ABSTRACT_JSON_OBJECT_ELEMENT_PARSER(json_int64(i));
  }
  bool json_uint64(uint64_t i) override
  {
    return ABSTRACT_JSON_OBJECT_ELEMENT_PARSER(json_uint64(i));
  }
  bool json_double(double d) override
  {
    return ABSTRACT_JSON_OBJECT_ELEMENT_PARSER(json_double(d));
  }
  bool json_raw_number(const Ch *str, SizeType length, bool copy) override
  {
    return ABSTRACT_JSON_OBJECT_ELEMENT_PARSER(RawNumber(str, length, copy));
  }
  bool json_string(const Ch *str, SizeType length, bool copy) override
  {
    return ABSTRACT_JSON_OBJECT_ELEMENT_PARSER(json_string(str, length, copy));
  }
  bool json_start_object() override
  {
    ++object_level;
    switch(state)
      {
      case Start:
        // TODO we also set skip in element_parsers constructors.
        // This makes interface more verbose.
        // Maybe we can set skip=true and (on_skipped=[]{}) by default,
        // if we always reset them on object/array start.
        // NB: this is somewhat less safe.
        reset_element_parsers(this->skip);
        state = Inside;
        return true;
      case Inside:
        return ABSTRACT_JSON_OBJECT_ELEMENT_PARSER(json_start_object());
      default: return base_type::json_start_object();
      }
  }
  bool json_key(const Ch *str, SizeType length, bool copy) override
  {
    switch(state)
      {
        case Inside: {
          if(object_level == 1 && array_level == 0)
            {
              key = std::string(str, length);
              return true;
            }
          return ABSTRACT_JSON_OBJECT_ELEMENT_PARSER(
            json_key(str, length, copy));
        }
      default: return base_type::json_key(str, length, copy);
      }
  }
  bool json_end_object(SizeType memberCount) override
  {
    --object_level;
    switch(state)
      {
        case Inside: {
          if(object_level == 0)
            {
              state = End;
              this->on_end();
              return true;
            }
          const bool element_res = ABSTRACT_JSON_OBJECT_ELEMENT_PARSER(
            json_end_object(memberCount));
          if(object_level == 1 && array_level == 0)
            key.clear();
          return element_res;
        }
      default: return base_type::json_end_object(memberCount);
      }
  }
  bool json_start_array() override
  {
    ++array_level;
    return ABSTRACT_JSON_OBJECT_ELEMENT_PARSER(json_start_array());
  }
  bool json_end_array(SizeType elementCount) override
  {
    --array_level;
    switch(state)
      {
        case Inside: {
          if(object_level == 0)
            {
              state = End;
              this->on_end();
              return true;
            }
          const bool element_res = ABSTRACT_JSON_OBJECT_ELEMENT_PARSER(
            json_end_array(elementCount));
          if(object_level == 1 && array_level == 0)
            key.clear();
          return element_res;
        }
      default: return base_type::json_end_array(elementCount);
      }
  }
};
template <class TResult>
Abstract_Json_Object_Parser<TResult>::Abstract_Json_Object_Parser(
  bool skip, const std::function<void(TResult &&)> &on_parsed,
  const std::function<void()> &on_skipped)
    : base_type(skip, on_parsed, on_skipped)
{}
template <class TResult>
void Abstract_Json_Object_Parser<TResult>::reset(const bool skip)
{
  this->skip = skip;
  this->state = Start;
  this->key.clear();
  this->object_level = 0;
  this->array_level = 0;
  reset_element_parsers(skip);
  clear_result();
}
