#include "../Damped_Rational_State.hxx"

void Damped_Rational_State::json_string(const std::string &s)
{
  if(parsing_constant)
    {
      constant_state.json_string(s);
      parsing_constant=false;
      value.constant=constant_state.value;
    }
  else if(parsing_base)
    {
      base_state.json_string(s);
      parsing_base=false;
      value.base=base_state.value;
    }
  else if(parsing_poles)
    {
      poles_state.json_string(s);
    }
  else
    {
      throw std::runtime_error("Invalid input file.  Unexpected string '" + s
                               + "' inside '" + name + "'");
    }
}
