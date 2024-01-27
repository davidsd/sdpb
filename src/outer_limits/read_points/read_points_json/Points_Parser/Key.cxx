#include "../Points_Parser.hxx"

bool Points_Parser::Key(const Ch *str, rapidjson::SizeType length, bool)
{
  std::string key(str, length);
  if(inside)
    {
      if(parsing_points)
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key
                                   , "' inside '", points_state.name , "'.");
        }
      else if(key == points_state.name)
        {
          parsing_points = true;
        }
      else
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key
                                   , "' inside the main object.");
        }
    }
  else
    {
      RUNTIME_ERROR("Found a key outside of the main object: "
                               + key);
    }
  return true;
}
