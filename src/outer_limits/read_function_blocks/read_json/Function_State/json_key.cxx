#include "../Function_State.hxx"

void Function_State::json_key(const std::string &key)
{
  if(parsing_max_delta)
    {
      RUNTIME_ERROR("Invalid input file. Unexpected key '", key, "' inside '",
                    name, "." + max_delta_state.name, "'.");
    }
  else if(parsing_epsilon_value)
    {
      RUNTIME_ERROR("Invalid input file. Unexpected key '", key, "' inside '",
                    name, "." + epsilon_value_state.name, "'.");
    }
  else if(parsing_infinity_value)
    {
      RUNTIME_ERROR("Invalid input file. Unexpected key '", key, "' inside '",
                    name, "." + infinity_value_state.name, "'.");
    }
  else if(parsing_chebyshev_values)
    {
      RUNTIME_ERROR("Invalid input file. Unexpected key '", key, "' inside '",
                    name, "." + chebyshev_values_state.name, "'.");
    }
  else if(key == max_delta_state.name)
    {
      parsing_max_delta = true;
    }
  else if(key == epsilon_value_state.name)
    {
      parsing_epsilon_value = true;
    }
  else if(key == infinity_value_state.name)
    {
      parsing_infinity_value = true;
    }
  else if(key == chebyshev_values_state.name)
    {
      parsing_chebyshev_values = true;
    }
  else
    {
      RUNTIME_ERROR("Invalid input file. Unexpected key '", key, "' inside '",
                    name, "'");
    }
}
