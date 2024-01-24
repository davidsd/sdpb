#include "../Checkpoint_Parser.hxx"
#include "sdpb_util/assert.hxx"

bool Checkpoint_Parser::EndObject(rapidjson::SizeType)
{
  if(inside)
    {
      if(parsing_generation)
        {
          RUNTIME_ERROR(
            "Invalid input file. Unexpected object ending inside '",
            generation_state.name, "'");
        }
      else if(parsing_threshold)
        {
          RUNTIME_ERROR(
            "Invalid input file. Unexpected object ending inside '",
            threshold_state.name, "'");
        }
      else if(parsing_c_scale)
        {
          RUNTIME_ERROR(
            "Invalid input file. Unexpected object ending inside '",
            c_scale_state.name, "'");
        }
      else if(parsing_yp)
        {
          RUNTIME_ERROR(
            "Invalid input file. Unexpected object ending inside '",
            yp_state.name, "'");
        }
      else if(parsing_b)
        {
          RUNTIME_ERROR(
            "Invalid input file. Unexpected object ending inside '",
            b_state.name, "'");
        }
      else if(parsing_y_transform)
        {
          RUNTIME_ERROR(
            "Invalid input file. Unexpected object ending inside '",
            y_transform_state.name, "'");
        }
      else if(parsing_points)
        {
          RUNTIME_ERROR(
            "Invalid input file. Unexpected object ending inside '",
            points_state.name, "'");
        }
      else
        {
          inside = false;
        }
    }
  else
    {
      RUNTIME_ERROR("Invalid input file. "
                    "Unexpected object ending outside of the main object");
    }
  return true;
}
