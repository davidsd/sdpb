#pragma once
#include <functional>

template <class TValue> struct Abstract_Element_Parser
{
  using value_type = TValue;

protected:
  bool skip;

public:
  void on_end()
  {
    if(this->skip)
      on_skipped_func();
    else
      on_parsed_func(get_result());
  }

  virtual value_type get_result() = 0;
  virtual void reset(bool skip) = 0;

protected:
  std::function<void(value_type &&)> on_parsed_func;
  std::function<void()> on_skipped_func;

protected:
  // TODO reorder and add default values: on_parsed, skip=false, on_skipped=[]{}
  // NB: if we set skip=false by default, then users (e.g. Array_Parser)
  // are responsible to set correct value before starting.
  Abstract_Element_Parser(const bool skip,
                          const std::function<void(value_type &&)> &on_parsed,
                          const std::function<void()> &on_skipped)
      : skip(skip), on_parsed_func(on_parsed), on_skipped_func(on_skipped)
  {}
  ~Abstract_Element_Parser() = default;
};
