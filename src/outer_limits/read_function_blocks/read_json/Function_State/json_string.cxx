#include "../Function_State.hxx"

void Function_State::json_string(const std::string &s)
{
  if(parsing_max_delta)
    {
      max_delta_state.json_string(s);
      max_delta_mpfr = Boost_Float(s);
      parsing_max_delta = false;
      value.max_delta = max_delta_state.value;
    }
  else if(parsing_epsilon_value)
    {
      epsilon_value_state.json_string(s);
      parsing_epsilon_value = false;
      value.epsilon_value = epsilon_value_state.value;
    }
  else if(parsing_infinity_value)
    {
      infinity_value_state.json_string(s);
      parsing_infinity_value = false;
      value.infinity_value = infinity_value_state.value;
    }
  else if(parsing_chebyshev_values)
    {
      chebyshev_values_state.json_string(s);
    }
  else
    {
      RUNTIME_ERROR("Invalid input file. Unexpected string '" , s
                               , "' inside '" , name , "'");
    }
}
