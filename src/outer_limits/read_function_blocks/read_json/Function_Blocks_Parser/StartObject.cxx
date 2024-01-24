#include "../Function_Blocks_Parser.hxx"
#include "sdpb_util/assert.hxx"

bool Function_Blocks_Parser::StartObject()
{
  if(inside)
    {
      if(parsing_objective)
        {
          RUNTIME_ERROR(
            "Invalid input file. Found an object inside '"
            , objective_state.name , "'");
        }
      else if(parsing_normalization)
        {
          RUNTIME_ERROR(
            "Invalid input file. Found an object inside '"
            , normalization_state.name , "'");
        }
      else if(parsing_functions)
        {
          functions_state.json_start_object();
        }
      else
        {
          RUNTIME_ERROR(
            "Invalid input file. Unknown object inside the main object");
        }
    }
  else
    {
      inside = true;
    }
  return true;
}
