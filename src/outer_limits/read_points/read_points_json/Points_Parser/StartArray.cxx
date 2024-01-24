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
          RUNTIME_ERROR(
            "Invalid input file. Unknown array inside the main object");
        }
    }
  else
    {
      RUNTIME_ERROR("Found an array outside of the SDP.");
    }
  return true;
}
