#include "../Checkpoint_Parser.hxx"
#include "sdpb_util/assert.hxx"

bool Checkpoint_Parser::Key(const Ch *str, rapidjson::SizeType length, bool)
{
  std::string key(str, length);
  if(inside)
    {
      if(parsing_generation)
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key,
                        "' inside '", generation_state.name, "'.");
        }
      else if(parsing_threshold)
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key,
                        "' inside '", threshold_state.name, "'.");
        }
      else if(parsing_c_scale)
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key,
                        "' inside '", c_scale_state.name , "'.");
        }
      else if(parsing_yp)
        {
          RUNTIME_ERROR("Invalid input file. Found the key '" , key
                                   , "' inside '" , yp_state.name , "'.");
        }
      else if(parsing_b)
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key
                                   , "' inside '" , b_state.name , "'.");
        }
      else if(parsing_y_transform)
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key
                                   , "' inside '" , y_transform_state.name
                                   , "'.");
        }
      else if(parsing_points)
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key
                                   , "' inside '" , points_state.name , "'.");
        }

      else if(key == generation_state.name)
        {
          parsing_generation = true;
        }
      else if(key == threshold_state.name)
        {
          parsing_threshold = true;
        }
      else if(key == c_scale_state.name)
        {
          parsing_c_scale = true;
        }
      else if(key == yp_state.name)
        {
          parsing_yp = true;
        }
      else if(key == b_state.name)
        {
          parsing_b = true;
        }
      else if(key == y_transform_state.name)
        {
          parsing_y_transform = true;
        }
      else if(key == points_state.name)
        {
          parsing_points = true;
        }
      else
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key
                                   + "' inside the main object.");
        }
    }
  else
    {
      RUNTIME_ERROR("Found a key outside of the main object: "
                               + key);
    }
  return true;
}
