#include "../Checkpoint_Parser.hxx"
#include "sdpb_util/assert.hxx"

bool Checkpoint_Parser::StartObject()
{
  if(inside)
    {
      if(parsing_generation)
        {
          RUNTIME_ERROR("Invalid input file. Unexpected object start inside '",
                        generation_state.name, "'");
        }
      else if(parsing_threshold)
        {
          RUNTIME_ERROR("Invalid input file. Unexpected object start inside '",
                        threshold_state.name, "'");
        }
      else if(parsing_c_scale)
        {
          RUNTIME_ERROR("Invalid input file. Unexpected object start inside '",
                        c_scale_state.name, "'");
        }
      else if(parsing_yp)
        {
          RUNTIME_ERROR("Invalid input file. Unexpected object start inside '",
                        yp_state.name, "'");
        }
      else if(parsing_b)
        {
          RUNTIME_ERROR("Invalid input file. Unexpected object start inside '",
                        b_state.name, "'");
        }
      else if(parsing_y_transform)
        {
          RUNTIME_ERROR("Invalid input file. Unexpected object start inside '",
                        y_transform_state.name, "'");
        }
      else if(parsing_points)
        {
          RUNTIME_ERROR("Invalid input file. Unexpected object start inside '",
                        points_state.name, "'");
        }
      else
        {
          RUNTIME_ERROR(
            "Invalid input file. Unknown object inside the main object");
        }
    }
  else
    {
      inside = true;
    }
  return true;
}
