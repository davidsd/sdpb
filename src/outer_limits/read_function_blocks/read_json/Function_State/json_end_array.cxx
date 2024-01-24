#include "../Function_State.hxx"
#include "sdpb_util/assert.hxx"

void Function_State::json_end_array()
{
  if(parsing_max_delta)
    {
      RUNTIME_ERROR("Invalid input file. Unexpected array end inside '", name,
                    ".", max_delta_state.name, "'.");
    }
  else if(parsing_epsilon_value)
    {
      RUNTIME_ERROR("Invalid input file. Unexpected array end inside '", name,
                    ".", epsilon_value_state.name, "'.");
    }
  else if(parsing_infinity_value)
    {
      RUNTIME_ERROR("Invalid input file. Unexpected array end inside '", name,
                    ".", infinity_value_state.name, "'.");
    }
  else if(parsing_chebyshev_values)
    {
      chebyshev_values_state.json_end_array();
      parsing_chebyshev_values = chebyshev_values_state.inside;
    }
  else
    {
      RUNTIME_ERROR("Invalid input file. Unexpected array end inside '", name,
                    "'");
    }
}
