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
//
// In case of a single polynomial, one can omit "product"
// and use Polynomial_Power directly, e.g.:
// {
//   "polynomial": ["0.2", "1.0"],
//   "power": "0.3"
// }
class Json_Polynomial_Power_Product_Parser final
    : public Abstract_Json_Object_Parser<Polynomial_Power_Product>
{
private:
  // terms entering the product
  std::optional<std::vector<Polynomial_Power>> terms;
  // special case: single Polynomial_Power
  std::optional<std::vector<Boost_Float>> single_polynomial;
  std::optional<Boost_Float> single_power;

  Json_Vector_Parser<Json_Polynomial_Power_Parser> polynomial_powers_parser;
  Json_Float_Vector_Parser<Boost_Float> polynomial_parser;
  Json_Float_Parser<Boost_Float> power_parser;

public:
  Json_Polynomial_Power_Product_Parser(
    const bool skip,
    const std::function<void(Polynomial_Power_Product &&)> &on_parsed,
    const std::function<void()> &on_skipped = [] {})
      : Abstract_Json_Object_Parser(skip, on_parsed, on_skipped),
        polynomial_powers_parser(
          skip,
          [this](std::vector<Polynomial_Power> &&value) {
            this->terms = std::move(value);
          }),
        polynomial_parser(skip,
                          [this](std::vector<Boost_Float> &&value) {
                            this->single_polynomial = std::move(value);
                          }),
        power_parser(skip, [this](Boost_Float &&value) {
          this->single_power = std::move(value);
        })
  {}

protected:
  Abstract_Json_Reader_Handler &element_parser(const std::string &key) override
  {
    if(key == "product")
      return polynomial_powers_parser;
    if(key == "polynomial")
      return polynomial_parser;
    if(key == "power")
      return power_parser;
    RUNTIME_ERROR("Json_Polynomial_Power_Product_Parser: unknown key=", key);
  }

public:
  value_type get_result() override
  {
    Polynomial_Power_Product result;
    if(terms.has_value())
      {
        ASSERT(!single_polynomial.has_value(),
               "Json_Polynomial_Power_Product_Parser: Cannot specify both "
               "\"product\" and \"polynomial\"");
        ASSERT(!single_power.has_value(),
               "Json_Polynomial_Power_Product_Parser: Cannot specify both "
               "\"product\" and \"power\"");
        result.terms = std::move(terms.value());
      }
    else
      {
        ASSERT(single_polynomial.has_value() && single_power.has_value(),
               "Json_Polynomial_Power_Product_Parser: expected either "
               "a \"product\" or a single \"polynomial\" raised to \"power\"");
        auto& pp = result.terms.emplace_back();
        pp.polynomial = std::move(single_polynomial.value());
        pp.power = std::move(single_power.value());
      }
    return result;
  }
  void clear_result() override
  {
    terms.reset();
    single_polynomial.reset();
    single_power.reset();
  }
  void reset_element_parsers(const bool skip) override
  {
    polynomial_powers_parser.reset(skip);
    polynomial_parser.reset(skip);
    power_parser.reset(skip);
  }
};
