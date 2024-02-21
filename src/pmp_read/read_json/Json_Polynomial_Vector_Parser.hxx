#pragma once

#include "Json_Polynomial_Parser.hxx"
#include "sdpb_util/json/Json_Float_Parser.hxx"
#include "sdpb_util/json/Json_Vector_Parser.hxx"
#include "pmp/Polynomial.hxx"

class Json_Polynomial_Vector_Parser final
    : public Abstract_Json_Vector_Parser<Polynomial_Vector,
                                         Json_Polynomial_Parser>
{
  Polynomial_Vector result;

public:
  using element_type = Polynomial;
  using value_type = Polynomial_Vector;

  Json_Polynomial_Vector_Parser(
    const bool skip,
    const std::function<void(Polynomial_Vector &&)> &on_parsed,
    const std::function<void()> &on_skipped = [] {})
      : Abstract_Json_Vector_Parser(skip, on_parsed, on_skipped)
  {}

  void clear_result() override { result.clear(); }
  void on_element_parsed(element_type &&value, size_t index) override
  {
    ASSERT_EQUAL(index, result.size());
    result.push_back(value);
  }
  void on_element_skipped(size_t /*index*/) override {}
  value_type get_result() override { return std::move(result); }
};
