#include "../Function_Blocks_Parser.hxx"

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
          throw std::runtime_error(
            "Invalid input file.  Unknown array inside the main object");
        }
    }
  else
    {
      throw std::runtime_error("Found an array outside of the SDP.");
    }
  return true;
}
