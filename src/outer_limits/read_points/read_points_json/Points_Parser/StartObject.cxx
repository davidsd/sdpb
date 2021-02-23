#include "../Points_Parser.hxx"

bool Points_Parser::StartObject()
{
  if(inside)
    {
      if(parsing_points)
        {
          throw std::runtime_error(
            "Invalid input file.  Found an object inside '" + points_state.name
            + "'");
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
