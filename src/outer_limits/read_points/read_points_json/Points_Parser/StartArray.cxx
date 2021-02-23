#include "../Points_Parser.hxx"

bool Points_Parser::StartArray()
{
  if(inside)
    {
      if(parsing_points)
        {
          points_state.json_start_array();
        }
      else
        {
          throw std::runtime_error(
            "Invalid input file.  Unknown array inside the main object");
        }
    }
  else
    {
      throw std::runtime_error("Found an array outside of the SDP.");
    }
  return true;
}
