#include "../Checkpoint_Parser.hxx"

bool Checkpoint_Parser::StartObject()
{
  if(inside)
    {
      if(parsing_generation)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected object start inside '"
            + generation_state.name + "'");
        }
      else if(parsing_threshold)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected object start inside '"
            + threshold_state.name + "'");
        }
      else if(parsing_c_scale)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected object start inside '"
            + c_scale_state.name + "'");
        }
      else if(parsing_yp)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected object start inside '"
            + yp_state.name + "'");
        }
      else if(parsing_b)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected object start inside '"
            + b_state.name + "'");
        }
      else if(parsing_y_transform)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected object start inside '"
            + y_transform_state.name + "'");
        }
      else if(parsing_points)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected object start inside '"
            + points_state.name + "'");
        }
      else
        {
          throw std::runtime_error(
            "Invalid input file.  Unknown object inside the main object");
        }
    }
  else
    {
      inside = true;
    }
  return true;
}
