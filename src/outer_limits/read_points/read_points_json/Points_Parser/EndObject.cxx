#include "../Points_Parser.hxx"

bool Points_Parser::EndObject(rapidjson::SizeType)
{
  if(inside)
    {
      if(parsing_points)
        {
          RUNTIME_ERROR(
            "Invalid input file. Found an object end inside '"
            + points_state.name , "'");
        }
      else
        {
          inside = false;
        }
    }
  else
    {
      RUNTIME_ERROR(
        "Invalid input file. Ending an object outside of the SDP");
    }
  return true;
}
