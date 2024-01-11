#include "../JSON_Parser.hxx"

bool JSON_Parser::StartObject()
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
      else if(parsing_positive_matrices_with_prefactor)
        {
          positive_matrices_with_prefactor_state.json_start_object();
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
