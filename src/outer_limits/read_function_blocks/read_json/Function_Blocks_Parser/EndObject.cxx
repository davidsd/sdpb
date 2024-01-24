#include "../Function_Blocks_Parser.hxx"
#include "sdpb_util/assert.hxx"

bool Function_Blocks_Parser::EndObject(rapidjson::SizeType)
{
  if(inside)
    {
      if(parsing_objective)
        {
          RUNTIME_ERROR("Invalid input file. Found an object end inside '",
                        objective_state.name, "'");
        }
      else if(parsing_normalization)
        {
          RUNTIME_ERROR("Invalid input file. Found an object end inside '",
                        normalization_state.name, "'");
        }
      else if(parsing_functions)
        {
          functions_state.json_end_object();
          parsing_functions = functions_state.inside;
        }
      else
        {
          inside = false;
        }
    }
  else
    {
      RUNTIME_ERROR("Invalid input file. Ending an object outside of the SDP");
    }
  return true;
}
