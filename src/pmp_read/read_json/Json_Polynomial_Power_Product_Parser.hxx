#pragma once

#include "Json_Polynomial_Power_Parser.hxx"
#include "pmp/Polynomial_Power_Product.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/json/Abstract_Json_Object_Parser.hxx"
#include "sdpb_util/json/Json_Float_Parser.hxx"

// Product of polynomials raised to arbitrary powers,
// e.g. (0.2 + x)^0.3 * (10.1 + x + 3.1 * x^2)^0.8 :
// {
//   "product" : [
//     {"polynomial": ["0.2", "1.0"], "power": "0.3"},
//     {"polynomial": ["10.1", "1.0", "3.1"], "power": "0.8"}
//   ]
// }
class Json_Polynomial_Power_Product_Parser final
    : public Abstract_Json_Object_Parser<Polynomial_Power_Product>
{
private:
  Polynomial_Power_Product result;
  Json_Vector_Parser<Json_Polynomial_Power_Parser> polynomial_powers_parser;

public:
  Json_Polynomial_Power_Product_Parser(
    const bool skip,
    const std::function<void(Polynomial_Power_Product &&)> &on_parsed,
    const std::function<void()> &on_skipped = [] {})
      : Abstract_Json_Object_Parser(skip, on_parsed, on_skipped),
        polynomial_powers_parser(
          skip, [this](std::vector<Polynomial_Power> &&value) {
            this->result.terms = std::move(value);
          })
  {}

protected:
  Abstract_Json_Reader_Handler &element_parser(const std::string &key) override
  {
    if(key == "product")
      return polynomial_powers_parser;
    RUNTIME_ERROR("Json_Polynomial_Power_Product_Parser: unknown key=", key);
  }

public:
  value_type get_result() override { return std::move(result); }
  void clear_result() override { result.terms.clear(); }
  void reset_element_parsers(const bool skip) override
  {
    polynomial_powers_parser.reset(skip);
  }
};
