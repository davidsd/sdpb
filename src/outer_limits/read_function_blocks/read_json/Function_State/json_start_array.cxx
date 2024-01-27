#include "../Function_State.hxx"

void Function_State::json_start_array()
{
  if(parsing_max_delta)
    {
      RUNTIME_ERROR("Invalid input file. Unexpected array inside '"
                               + name , "." , max_delta_state.name , "'.");
    }
  else if(parsing_epsilon_value)
    {
      RUNTIME_ERROR("Invalid input file. Unexpected array inside '"
                               + name , "." , epsilon_value_state.name , "'.");
    }
  else if(parsing_infinity_value)
    {
      RUNTIME_ERROR("Invalid input file. Unexpected array inside '"
                               + name , "." , infinity_value_state.name
                               , "'.");
    }
  else if(parsing_chebyshev_values)
    {
      chebyshev_values_state.json_start_array();
    }
  else
    {
      RUNTIME_ERROR("Invalid input file. Unexpected array inside '"
                               + name , "'");
    }
}
