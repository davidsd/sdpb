#include "../Checkpoint_Parser.hxx"

bool Checkpoint_Parser::String(const Ch *str, rapidjson::SizeType length, bool)
{
  std::string s(str, length);
  if(inside)
    {
      if(parsing_generation)
        {
          generation_state.json_string(s);
          parsing_generation = false;
        }
      else if(parsing_threshold)
        {
          threshold_state.json_string(s);
          parsing_threshold = false;
        }
      else if(parsing_c_scale)
        {
          c_scale_state.json_string(s);
          parsing_c_scale = false;
        }
      else if(parsing_y)
        {
          y_state.json_string(s);
        }
      else if(parsing_b)
        {
          b_state.json_string(s);
        }
      else if(parsing_y_transform)
        {
          y_transform_state.json_string(s);
        }
      else if(parsing_points)
        {
          points_state.json_string(s);
        }
      else
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected string in the main object: '" + s
            + "'");
        }
    }
  else
    {
      throw std::runtime_error("Found a string outside of the SDP: '" + s
                               + "'");
    }
  return true;
}
