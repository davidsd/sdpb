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
  else if(parsing_polynomial)
    {
      polynomial_state.json_end_array();
      parsing_polynomial = polynomial_state.inside;
    }
  else
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected array end inside '" + name + "'");
    }
}
