#include "Block_Parser.hxx"
#include "sdpb_util/assert.hxx"

bool Block_Parser::Key(const Ch *str, rapidjson::SizeType length, bool)
{
  std::string key(str, length);
  if(inside)
    {
      if(parsing_bilinear_bases_even)
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key,
                        "' inside '", bilinear_bases_even_state.name, "'.");
        }
      else if(parsing_bilinear_bases_odd)
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key,
                        "' inside '", bilinear_bases_odd_state.name, "'.");
        }
      else if(parsing_c)
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key,
                        "' inside '", c_state.name, "'.");
        }
      else if(parsing_B)
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key,
                        "' inside '", B_state.name, "'.");
        }
      if(key == bilinear_bases_even_state.name)
        {
          parsing_bilinear_bases_even = true;
        }
      else if(key == bilinear_bases_odd_state.name)
        {
          parsing_bilinear_bases_odd = true;
        }
      else if(key == c_state.name)
        {
          parsing_c = true;
        }
      else if(key == B_state.name)
        {
          parsing_B = true;
        }
      else
        {
          RUNTIME_ERROR("Invalid input file. Found the key '", key,
                        "' inside the main object.");
        }
    }
  else
    {
      RUNTIME_ERROR("Found a key outside of the main object: ", key);
    }
  return true;
}
