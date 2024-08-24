#pragma once

#include "pmp/Polynomial_Power_Product.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/json/Abstract_Json_Object_Parser.hxx"
#include "sdpb_util/json/Json_Float_Parser.hxx"

// Polynomial raised to an arbitrary power,
// e.g. (0.2 + x)^0.3 :
// {
//   "polynomial" : ["0.2", "1.0"], "power" : "0.3"
// }

class Json_Polynomial_Power_Parser final
    : public Abstract_Json_Object_Parser<Polynomial_Power>
{
private:
  Polynomial_Power result;
  Json_Float_Vector_Parser<Boost_Float> polynomial_parser;
  Json_Float_Parser<Boost_Float> power_parser;

public:
  Json_Polynomial_Power_Parser(
    const bool skip, const std::function<void(Polynomial_Power &&)> &on_parsed,
    const std::function<void()> &on_skipped = [] {})
      : Abstract_Json_Object_Parser(skip, on_parsed, on_skipped),
        polynomial_parser(skip,
                    [this](std::vector<Boost_Float> &&value) {
                      this->result.polynomial = std::move(value);
                    }),
        power_parser(skip,
                        [this](Boost_Float &&value) {
                          this->result.power = std::move(value);
                        })
  {}

protected:
  Abstract_Json_Reader_Handler &element_parser(const std::string &key) override
  {
    if(key == "polynomial")
      return polynomial_parser;
    if(key == "power")
      return power_parser;
    RUNTIME_ERROR("Json_Polynomial_Power_Parser: unknown key=", key);
  }

public:
  value_type get_result() override { return std::move(result); }
  void clear_result() override
  {
    result.polynomial = {};
    result.power = 1;
  }
  void reset_element_parsers(const bool skip) override
  {
    polynomial_parser.reset(skip);
    power_parser.reset(skip);
  }
};
