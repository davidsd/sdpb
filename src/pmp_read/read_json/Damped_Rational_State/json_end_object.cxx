#include "../Damped_Rational_State.hxx"

void Damped_Rational_State::json_end_object()
{
  if(parsing_constant)
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected object end inside '" + name + "."
        + constant_state.name + "'.");
    }
  else if(parsing_base)
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected object end inside '" + name + "."
        + base_state.name + "'.");
    }
  else if(parsing_poles)
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected object end inside '" + name + "."
        + poles_state.name + "'.");
    }
  inside = false;
}
