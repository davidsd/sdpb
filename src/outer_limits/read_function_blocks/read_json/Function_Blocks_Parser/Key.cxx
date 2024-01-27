#include "../Function_Blocks_Parser.hxx"
#include "sdpb_util/assert.hxx"

bool Function_Blocks_Parser::Key(const Ch *str, rapidjson::SizeType length,
                                 bool)
{
  const std::string key(str, length);
  if(inside)
    {
      if(parsing_objective)
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key,
                        "' inside '", objective_state.name, "'.");
        }
      else if(parsing_normalization)
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key,
                        "' inside '", normalization_state.name, "'.");
        }
      else if(parsing_functions)
        {
          functions_state.json_key(key);
        }
      else if(key == objective_state.name)
        {
          parsing_objective = true;
        }
      else if(key == normalization_state.name)
        {
          parsing_normalization = true;
        }
      else if(key == functions_state.name)
        {
          parsing_functions = true;
        }
      else
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key,
                        "' inside the main object.");
        }
    }
  else
    {
      RUNTIME_ERROR("Found a key outside of the main object: ", key);
    }
  return true;
}
