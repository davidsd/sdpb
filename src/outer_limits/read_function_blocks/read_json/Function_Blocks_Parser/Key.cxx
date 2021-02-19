#include "../Function_Blocks_Parser.hxx"

bool Function_Blocks_Parser::Key(const Ch *str, rapidjson::SizeType length,
                                 bool)
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
      else if(parsing_points)
        {
          throw std::runtime_error("Invalid input file.  Found the key '" + key
                                   + "' inside '" + points_state.name
                                   + "'.");
        }
      else if(parsing_functions)
        {
          throw std::runtime_error("Invalid input file.  Found the key '" + key
                                   + "' inside '" + functions_state.name
                                   + "'.");
        }
      else if(key == objective_state.name)
        {
          parsing_objective = true;
        }
      else if(key == normalization_state.name)
        {
          parsing_normalization = true;
        }
      else if(key == points_state.name)
        {
          parsing_points = true;
        }
      else if(key == functions_state.name)
        {
          parsing_functions = true;
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
