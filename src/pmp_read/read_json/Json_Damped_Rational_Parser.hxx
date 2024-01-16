#pragma once

#include "pmp/Damped_Rational.hxx"
#include "sdpb_util/json/Abstract_Json_Object_Parser.hxx"
#include "sdpb_util/json/Json_Float_Parser.hxx"

// "DampedRational":
// {
//   "type": "object",
//   "properties":
//   {
//     "base": { "type": "string" },
//     "constant": { "type": "string" },
//     "poles":
//     {
//       "type": "array",
//       "items": { "type": "string" }
//     }
//   }
// },
class Json_Damped_Rational_Parser final
    : public Abstract_Json_Object_Parser<Damped_Rational>
{
private:
  Damped_Rational result;
  Json_Float_Parser<Boost_Float> base_parser;
  Json_Float_Parser<Boost_Float> constant_parser;
  Json_Float_Vector_Parser<Boost_Float> poles_parser;

public:
  Json_Damped_Rational_Parser(
    const bool skip, const std::function<void(Damped_Rational &&)> &on_parsed,
    const std::function<void()> &on_skipped = [] {})
      : Abstract_Json_Object_Parser(skip, on_parsed, on_skipped),
        // TODO set on_skipped = []{} everywhere by default
        base_parser(skip,
                    [this](Boost_Float &&value) {
                      this->result.base = std::move(value);
                    }),
        constant_parser(skip,
                        [this](Boost_Float &&value) {
                          this->result.constant = std::move(value);
                        }),
        poles_parser(skip, [this](std::vector<Boost_Float> &&value) {
          this->result.poles = std::move(value);
        })
  {}

protected:
  Abstract_Json_Reader_Handler &element_parser(const std::string &key) override
  {
    if(key == "base")
      return base_parser;
    if(key == "constant")
      return constant_parser;
    if(key == "poles")
      return poles_parser;
    El::RuntimeError("Json_Damped_Rational_Parser: unknown key=", key);
  }

public:
  value_type get_result() override { return std::move(result); }
  void clear_result() override
  {
    result.base.assign(0);
    result.constant.assign(0);
    result.poles.clear();
  }
  void reset_element_parsers(const bool skip) override
  {
    base_parser.reset(skip);
    constant_parser.reset(skip);
    poles_parser.reset(skip);
  }
};
