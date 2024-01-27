#pragma once
#include "Abstract_Element_Parser.hxx"
#include "Abstract_Json_Reader_Handler.hxx"

#include <functional>

template <class TValue>
class Abstract_Json_Element_Parser
    : public Abstract_Element_Parser<TValue>,
      public Abstract_Json_Reader_Handler
{
public:
  using value_type = typename Abstract_Element_Parser<TValue>::value_type;

protected:
  Abstract_Json_Element_Parser(
    const bool skip, const std::function<void(value_type &&)> &on_parsed,
    const std::function<void()> &on_skipped)
      : Abstract_Element_Parser<TValue>(skip, on_parsed, on_skipped)
  {}
  ~Abstract_Json_Element_Parser() = default;
};
