#include "../Points_Parser.hxx"

bool Points_Parser::EndObject(rapidjson::SizeType)
{
  if(inside)
    {
      if(parsing_points)
        {
          throw std::runtime_error(
            "Invalid input file.  Found an object end inside '"
            + points_state.name + "'");
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
