#include "../Damped_Rational_State.hxx"

void Damped_Rational_State::json_end_array()
{
  if(parsing_constant)
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected array end inside '" + name + "."
        + constant_state.name + "'.");
    }
  else if(parsing_base)
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected array end inside '" + name + "."
        + base_state.name + "'.");
    }
  else if(parsing_poles)
    {
      poles_state.json_end_array();
      parsing_poles = poles_state.inside;
      if(!poles_state.inside)
        {
          std::swap(value.poles,poles_state.value);
        }
    }
  else
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected array end inside '" + name + "'");
    }
}
