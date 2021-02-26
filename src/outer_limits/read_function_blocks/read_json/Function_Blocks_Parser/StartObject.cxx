#include "../Function_Blocks_Parser.hxx"

bool Function_Blocks_Parser::StartObject()
{
  if(inside)
    {
      if(parsing_objective)
        {
          throw std::runtime_error(
            "Invalid input file.  Found an object inside '"
            + objective_state.name + "'");
        }
      else if(parsing_normalization)
        {
          throw std::runtime_error(
            "Invalid input file.  Found an object inside '"
            + normalization_state.name + "'");
        }
      else if(parsing_functions)
        {
          functions_state.json_start_object();
        }
      else
        {
          throw std::runtime_error(
            "Invalid input file.  Unknown object inside the main object");
        }
    }
  else
    {
      inside = true;
    }
  return true;
}
