#include "../JSON_Parser.hxx"

bool JSON_Parser::EndObject(rapidjson::SizeType)
{
  if(inside)
    {
      if(parsing_objective)
        {
          throw std::runtime_error(
            "Invalid input file.  Found an object end inside '"
            + objective_state.name + "'");
        }
      else if(parsing_normalization)
        {
          throw std::runtime_error(
            "Invalid input file.  Found an object end inside '"
            + normalization_state.name + "'");
        }
      else if(parsing_positive_matrices_with_prefactor)
        {
          positive_matrices_with_prefactor_state.json_end_object();
        }
      else
        {
          inside = false;
        }
    }
  else
    {
      throw std::runtime_error(
        "Invalid input file.  Ending an object outside of the SDP");
    }
  return true;
}
