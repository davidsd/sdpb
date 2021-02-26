#include "../Function_State.hxx"

void Function_State::json_end_object()
{
  std::cout << "end\n" << std::flush;
  if(parsing_max_delta)
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected object end inside '" + name + "."
        + max_delta_state.name + "'.");
    }
  else if(parsing_infinity_value)
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected object end inside '" + name + "."
        + infinity_value_state.name + "'.");
    }
  else if(parsing_chebyshev_values)
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected object end inside '" + name + "."
        + chebyshev_values_state.name + "'.");
    }
  inside = false;
}
