#include "../Damped_Rational_State.hxx"

void Damped_Rational_State::json_start_array()
{
  if(parsing_constant)
    {
      throw std::runtime_error("Invalid input file.  Unexpected array inside '"
                               + name + "." + constant_state.name + "'.");
    }
  else if(parsing_base)
    {
      throw std::runtime_error("Invalid input file.  Unexpected array inside '"
                               + name + "." + base_state.name + "'.");
    }
  else if(parsing_poles)
    {
      poles_state.json_start_array();
    }
  else
    {
      throw std::runtime_error("Invalid input file.  Unexpected array inside '"
                               + name + "'");
    }
}
