#pragma once

#include "sdpb_util/json/Json_Float_Parser.hxx"
#include "sdpb_util/json/Json_Vector_Parser.hxx"
#include "pmp/Polynomial.hxx"

class Json_Polynomial_Parser final
    : public Abstract_Json_Vector_Parser<Polynomial,
                                         Json_Float_Parser<El::BigFloat>>
{
  Polynomial result{0, 0};

public:
  using element_type = El::BigFloat;
  using value_type = Polynomial;

  Json_Polynomial_Parser(
    const bool skip, const std::function<void(Polynomial &&)> &on_parsed,
    const std::function<void()> &on_skipped = [] {})
      : Abstract_Json_Vector_Parser(skip, on_parsed, on_skipped)
  {}

  void clear_result() override { result.coefficients.clear(); }
  void on_element_parsed(element_type &&value, size_t index) override
  {
    ASSERT(index == result.coefficients.size(), "index=", index,
           ", result.coefficients.size()=", result.coefficients.size());
    result.coefficients.push_back(std::forward<element_type>(value));
  }
  void on_element_skipped(size_t index) override {}
  value_type get_result() override { return std::move(result); }
};
