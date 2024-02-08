#pragma once

#include "outer_limits/Function.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/json/Abstract_Json_Object_Parser.hxx"
#include "sdpb_util/json/Json_Float_Parser.hxx"
#include "sdpb_util/Boost_Float.hxx"

using namespace std::string_literals;
struct Json_Function_Parser final : Abstract_Json_Object_Parser<Function>
{
private:
  Function result;
  Json_Float_Parser<El::BigFloat> max_delta_parser;
  Json_Float_Parser<El::BigFloat> epsilon_value_parser;
  Json_Float_Parser<El::BigFloat> infinity_value_parser;
  Json_Float_Vector_Parser<Boost_Float> chebyshev_values_parser;

  // Abstract_Json_Object_Parser implementation
public:
  Function get_result() override { return std::move(result); }

protected:
  Abstract_Json_Reader_Handler &element_parser(const std::string &key) override
  {
    if(key == "max_delta")
      return max_delta_parser;
    if(key == "epsilon_value")
      return epsilon_value_parser;
    if(key == "infinity_value")
      return infinity_value_parser;
    if(key == "chebyshev_values")
      return chebyshev_values_parser;
    RUNTIME_ERROR("Unknown key: '", key, "'");
  }

public:
  void reset_element_parsers(bool skip) override
  {
    max_delta_parser.reset(skip);
    epsilon_value_parser.reset(skip);
    infinity_value_parser.reset(skip);
    chebyshev_values_parser.reset(skip);
  }
  void clear_result() override
  {
    result.max_delta = 0;
    result.epsilon_value = 0;
    result.infinity_value = 0;
    result.chebyshev_coeffs.clear();
  }

public:
  Json_Function_Parser(const bool skip,
                       const std::function<void(Function &&)> &on_parsed,
                       const std::function<void()> &on_skipped)
      : Abstract_Json_Object_Parser(skip, on_parsed, on_skipped),
        max_delta_parser(skip,
                         [this](El::BigFloat &&max_delta) {
                           this->result.max_delta = std::move(max_delta);
                         }),
        epsilon_value_parser(skip,
                             [this](El::BigFloat &&epsilon_value) {
                               this->result.epsilon_value
                                 = std::move(epsilon_value);
                             }),
        infinity_value_parser(skip,
                              [this](El::BigFloat &&infinity_value) {
                                this->result.infinity_value
                                  = std::move(infinity_value);
                              }),
        chebyshev_values_parser(
          skip, [this](std::vector<Boost_Float> &&chebyshev_values) {
            // Convert from sampled values to chebyshev coefficients.
            const Boost_Float pi(boost::math::constants::pi<Boost_Float>());
            const size_t N(chebyshev_values.size());
            this->result.chebyshev_coeffs.resize(0);
            this->result.chebyshev_coeffs.reserve(N);
            std::vector<Boost_Float> coeffs(N, Boost_Float(0));
            for(size_t n(0); n < N; ++n)
              {
                Boost_Float coeff(0);
                for(size_t k(0); k < N; ++k)
                  {
                    coeff += 2
                             * cos((n * pi * (2 * (N - 1 - k) + 1)) / (2 * N))
                             * chebyshev_values[k] / N;
                  }
                this->result.chebyshev_coeffs.push_back(to_BigFloat(coeff));
              }
          })
  {}
};
