#include "../JSON_Parser.hxx"

bool JSON_Parser::EndArray(rapidjson::SizeType)
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
      else if(parsing_positive_matrices_with_prefactor)
        {
          positive_matrices_with_prefactor_state.json_end_array();
          parsing_positive_matrices_with_prefactor
            = positive_matrices_with_prefactor_state.inside;
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
