#include "../Points_Parser.hxx"

bool Points_Parser::String(const Ch *str, rapidjson::SizeType length, bool)
{
  std::string s(str, length);
  if(inside)
    {
      if(parsing_points)
        {
          points_state.json_string(s);
        }
      else
        {
          RUNTIME_ERROR(
            "Invalid input file. Unexpected string in the main object: '", s
            , "'");
        }
    }
  else
    {
      RUNTIME_ERROR("Found a string outside of the SDP: '", s
                               , "'");
    }
  return true;
}
