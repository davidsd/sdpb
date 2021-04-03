#include "../Checkpoint_Parser.hxx"

bool Checkpoint_Parser::EndObject(rapidjson::SizeType)
{
  if(inside)
    {
      if(parsing_generation)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected object ending inside '"
            + generation_state.name + "'");
        }
      else if(parsing_threshold)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected object ending inside '"
            + threshold_state.name + "'");
        }
      else if(parsing_c_scale)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected object ending inside '"
            + c_scale_state.name + "'");
        }
      else if(parsing_y)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected object ending inside '"
            + y_state.name + "'");
        }
      else if(parsing_b)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected object ending inside '"
            + b_state.name + "'");
        }
      else if(parsing_y_transform)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected object ending inside '"
            + y_transform_state.name + "'");
        }
      else if(parsing_points)
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected object ending inside '"
            + points_state.name + "'");
        }
      else
        {
          inside = false;
        }
    }
  else
    {
      throw std::runtime_error("Invalid input file.  Unexpected object ending "
                               "outside of the main object");
    }
  return true;
}
