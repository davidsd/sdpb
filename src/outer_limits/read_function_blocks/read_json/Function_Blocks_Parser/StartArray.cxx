#include "../Function_Blocks_Parser.hxx"
#include "sdpb_util/assert.hxx"

bool Function_Blocks_Parser::StartArray()
{
  if(inside)
    {
      if(parsing_objective)
        {
          objective_state.json_start_array();
        }
      else if(parsing_normalization)
        {
          normalization_state.json_start_array();
        }
      else if(parsing_functions)
        {
          functions_state.json_start_array();
        }
      else
        {
          RUNTIME_ERROR(
            "Invalid input file. Unknown array inside the main object");
        }
    }
  else
    {
      RUNTIME_ERROR("Found an array outside of the SDP.");
    }
  return true;
}
