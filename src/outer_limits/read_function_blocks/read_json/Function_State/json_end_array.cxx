#include "../Function_State.hxx"

void Function_State::json_end_array()
{
  if(parsing_max_delta)
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected array end inside '" + name + "."
        + max_delta_state.name + "'.");
    }
  else if(parsing_epsilon_value)
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected array end inside '" + name + "."
        + epsilon_value_state.name + "'.");
    }
  else if(parsing_infinity_value)
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected array end inside '" + name + "."
        + infinity_value_state.name + "'.");
    }
  else if(parsing_chebyshev_values)
    {
      chebyshev_values_state.json_end_array();
      parsing_chebyshev_values = chebyshev_values_state.inside;
    }
  else
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected array end inside '" + name + "'");
    }
}
