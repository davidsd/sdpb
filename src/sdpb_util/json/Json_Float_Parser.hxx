#pragma once

#include "Json_String_Element_Parser.hxx"
#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/assert.hxx"

#include <El.hpp>

#include <limits>

template <class TFloat>
class Json_Float_Parser final : public Json_String_Element_Parser<TFloat>
{
public:
  using typename Json_String_Element_Parser<TFloat>::value_type;

  Json_Float_Parser(
    const bool skip, const std::function<void(value_type &&)> &on_parsed,
    const std::function<void()> &on_skipped = [] {})
      : Json_String_Element_Parser<TFloat>(skip, on_parsed, on_skipped)
  {}

protected:
  TFloat from_string(const std::string &string_value) override
  {
    // GMP does not have inf or nan, so we approximate inf
    // with max double.
    // TODO throw error instead?
    using namespace std::string_literals;
    if(string_value == "inf"s)
      {
        constexpr auto inf = std::numeric_limits<double>::max();
        PRINT_WARNING("Replacing \"inf\" with ", inf);
        return TFloat(inf);
      }
    if(string_value == "-inf"s)
      {
        constexpr auto minus_inf = -std::numeric_limits<double>::max();
        PRINT_WARNING("Replacing \"-inf\" with ", minus_inf);
        return TFloat(minus_inf);
      }
    return TFloat(string_value);
  }
};

using Json_BigFloat_Parser = Json_Float_Parser<El::BigFloat>;
using Json_Boost_Float_Parser = Json_Float_Parser<Boost_Float>;
