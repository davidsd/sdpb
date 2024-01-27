#include "../Function_Blocks_Parser.hxx"
#include "sdpb_util/assert.hxx"

bool Function_Blocks_Parser::String(const Ch *str, rapidjson::SizeType length,
                                    bool)
{
  std::string s(str, length);
  if(inside)
    {
      if(parsing_objective)
        {
          objective_state.json_string(s);
        }
      else if(parsing_normalization)
        {
          normalization_state.json_string(s);
        }
      else if(parsing_functions)
        {
          functions_state.json_string(s);
        }
      else
        {
          RUNTIME_ERROR(
            "Invalid input file. Unexpected string in the main object: '" , s
            , "'");
        }
    }
  else
    {
      RUNTIME_ERROR("Found a string outside of the SDP: '" , s , "'");
    }
  return true;
}
