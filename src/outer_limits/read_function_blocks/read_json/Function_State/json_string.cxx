#include "../Function_State.hxx"

void Function_State::json_string(const std::string &s)
{
  if(parsing_max_delta)
    {
      max_delta_state.json_string(s);
      parsing_max_delta = false;
      value.max_delta = max_delta_state.value;
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
      throw std::runtime_error("Invalid input file.  Unexpected string '" + s
                               + "' inside '" + name + "'");
    }
}
