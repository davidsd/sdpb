#include "../Checkpoint_Parser.hxx"

bool Checkpoint_Parser::StartArray()
{
  if(inside)
    {
      if(parsing_generation)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected array start inside  '"
            + generation_state.name + "'");
        }
      else if(parsing_threshold)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected array start inside  '"
            + threshold_state.name + "'");
        }
      else if(parsing_c_scale)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected array start inside  '"
            + c_scale_state.name + "'");
        }
      else if(parsing_yp)
        {
          yp_state.json_start_array();
        }
      else if(parsing_b)
        {
          b_state.json_start_array();
        }
      else if(parsing_y_transform)
        {
          y_transform_state.json_start_array();
        }
      else if(parsing_points)
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
