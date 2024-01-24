#include "../Function_Blocks_Parser.hxx"
#include "sdpb_util/assert.hxx"

bool Function_Blocks_Parser::EndArray(rapidjson::SizeType)
{
  if(inside)
    {
      if(parsing_objective)
        {
          objective_state.json_end_array();
          parsing_objective = objective_state.inside;
        }
      else if(parsing_normalization)
        {
          normalization_state.json_end_array();
          parsing_normalization = normalization_state.inside;
        }
      else if(parsing_functions)
        {
          functions_state.json_end_array();
          parsing_functions = functions_state.inside;
        }
      else
        {
          RUNTIME_ERROR(
            "Invalid input file. Ending an array inside the main object.");
        }
    }
  else
    {
      RUNTIME_ERROR(
        "Invalid input file. Ending an array outside of the main object.");
    }
  return true;
}
