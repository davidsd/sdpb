#pragma once

#include "sdpb_util/Vector_State.hxx"
#include "sdpb_util/Number_State.hxx"
#include "sdp_read/Damped_Rational.hxx"

#include <vector>
#include <string>

using namespace std::string_literals;
struct Damped_Rational_State
{
  std::string name;
  bool inside = false, parsing_constant = false, parsing_base = false,
       parsing_poles = false;
  Damped_Rational value;
  Number_State<Boost_Float> constant_state, base_state;
  Vector_State<Number_State<Boost_Float>> poles_state;

  Damped_Rational_State(const std::vector<std::string> &names,
                        const size_t &offset)
      : name(names.at(offset)), constant_state("constant"s),
        base_state("base"s), poles_state({"poles"s, ""s})
  {}
  Damped_Rational_State(const std::initializer_list<std::string> &names)
      : Damped_Rational_State(names, 0)
  {}

  void json_key(const std::string &key);
  void json_string(const std::string &s);
  void json_start_array();
  void json_end_array();
  void json_start_object();
  void json_end_object();
};
