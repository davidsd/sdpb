#pragma once

#include "outer_limits/Function.hxx"
#include "sdpb_util/Vector_State.hxx"
#include "sdpb_util/Number_State.hxx"
#include "sdpb_util/Boost_Float.hxx"

using namespace std::string_literals;
struct Function_State
{
  std::string name;
  bool inside = false;
  bool parsing_max_delta = false, parsing_epsilon_value = false,
       parsing_infinity_value = false, parsing_chebyshev_values = false;
  Number_State<El::BigFloat> max_delta_state, epsilon_value_state,
    infinity_value_state;
  Vector_State<Number_State<Boost_Float>> chebyshev_values_state;

  Boost_Float max_delta_mpfr;
  Function value;

  Function_State(const std::vector<std::string> &names, const size_t &offset)
      : name(names.at(offset)), max_delta_state("max_delta"s),
        epsilon_value_state("epsilon_value"s),
        infinity_value_state("infinity_value"s),
        chebyshev_values_state({"chebyshev_values"s, ""s})
  {}

  void json_key(const std::string &key);
  void json_string(const std::string &s);
  void json_start_array();
  void json_end_array();
  void json_start_object();
  void json_end_object();
};
