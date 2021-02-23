#include "../Points_Parser.hxx"

bool Points_Parser::EndArray(rapidjson::SizeType)
{
  if(inside)
    {
      if(parsing_points)
        {
          points_state.json_end_array();
          parsing_points = points_state.inside;
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
