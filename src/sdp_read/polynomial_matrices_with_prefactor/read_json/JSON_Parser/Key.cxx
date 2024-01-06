#include "../JSON_Parser.hxx"

bool JSON_Parser::Key(const Ch *str, rapidjson::SizeType length, bool)
{
  std::string key(str, length);
  if(inside)
    {
      if(parsing_objective)
        {
          throw std::runtime_error("Invalid input file.  Found the key '" + key
                                   + "' inside '" + objective_state.name
                                   + "'.");
        }
      else if(parsing_normalization)
        {
          throw std::runtime_error("Invalid input file.  Found the key '" + key
                                   + "' inside '" + normalization_state.name
                                   + "'.");
        }
      else if(parsing_positive_matrices_with_prefactor)
        {
          positive_matrices_with_prefactor_state.json_key(key);
        }
      else if(key == objective_state.name)
        {
          parsing_objective = true;
        }
      else if(key == normalization_state.name)
        {
          parsing_normalization = true;
        }
      else if(key == positive_matrices_with_prefactor_state.name)
        {
          parsing_positive_matrices_with_prefactor = true;
        }
      else
        {
          throw std::runtime_error("Invalid input file.  Found the key '" + key
                                   + "' inside the main object.");
        }
    }
  else
    {
      throw std::runtime_error("Found a key outside of the main object: "
                               + key);
    }
  return true;
}
