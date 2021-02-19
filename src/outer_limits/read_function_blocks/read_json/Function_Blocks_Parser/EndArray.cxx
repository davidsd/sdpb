#include "../Function_Blocks_Parser.hxx"

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
      else if(parsing_points)
        {
          points_state.json_end_array();
          parsing_points = points_state.inside;
        }
      else if(parsing_functions)
        {
          functions_state.json_end_array();
          parsing_functions = functions_state.inside;
        }
      else
        {
          throw std::runtime_error(
            "Invalid input file.  Ending an array inside the main object.");
        }
    }
  else
    {
      throw std::runtime_error(
        "Invalid input file.  Ending an array outside of the main object.");
    }
  return true;
}
