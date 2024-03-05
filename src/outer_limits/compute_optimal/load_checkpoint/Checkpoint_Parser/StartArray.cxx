#include "../Checkpoint_Parser.hxx"
#include "sdpb_util/assert.hxx"

bool Checkpoint_Parser::StartArray()
{
  if(inside)
    {
      if(parsing_generation)
        {
          RUNTIME_ERROR(
            "Invalid input file. Unexpected array start inside  '",
            generation_state.name, "'");
        }
      else if(parsing_threshold)
        {
          RUNTIME_ERROR(
            "Invalid input file. Unexpected array start inside  '",
            threshold_state.name, "'");
        }
      else if(parsing_c_scale)
        {
          RUNTIME_ERROR(
            "Invalid input file. Unexpected array start inside  '",
            c_scale_state.name, "'");
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
