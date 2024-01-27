#include "../Checkpoint_Parser.hxx"
#include "sdpb_util/assert.hxx"

bool Checkpoint_Parser::EndArray(rapidjson::SizeType)
{
  if(inside)
    {
      if(parsing_generation)
        {
          RUNTIME_ERROR(
            "Invalid input file. Unexpected array ending inside  '",
            generation_state.name, "'");
        }
      else if(parsing_threshold)
        {
          RUNTIME_ERROR(
            "Invalid input file. Unexpected array ending inside  '",
            threshold_state.name, "'");
        }
      else if(parsing_c_scale)
        {
          RUNTIME_ERROR(
            "Invalid input file. Unexpected array ending inside  '",
            c_scale_state.name, "'");
        }
      else if(parsing_yp)
        {
          yp_state.json_end_array();
          parsing_yp = yp_state.inside;
        }
      else if(parsing_b)
        {
          b_state.json_end_array();
          parsing_b = b_state.inside;
        }
      else if(parsing_y_transform)
        {
          y_transform_state.json_end_array();
          parsing_y_transform = y_transform_state.inside;
        }
      else if(parsing_points)
        {
          points_state.json_end_array();
          parsing_points = points_state.inside;
        }
      else
        {
          RUNTIME_ERROR("Invalid input file. "
                        "Unexpected array ending inside the main object.");
        }
    }
  else
    {
      RUNTIME_ERROR("Invalid input file. "
                    " Unexpected array ending outside of the main object.");
    }
  return true;
}
