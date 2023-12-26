#include "Block_Parser.hxx"

bool Block_Parser::Key(const Ch *str, rapidjson::SizeType length, bool)
{
  std::string key(str, length);
  if(inside)
    {
      if(parsing_bilinear_bases_even)
        {
          throw std::runtime_error("Invalid input file.  Found the key '" + key
                                   + "' inside '"
                                   + bilinear_bases_even_state.name + "'.");
        }
      else if(parsing_bilinear_bases_odd)
        {
          throw std::runtime_error("Invalid input file.  Found the key '" + key
                                   + "' inside '"
                                   + bilinear_bases_odd_state.name + "'.");
        }
      else if(parsing_c)
        {
          throw std::runtime_error("Invalid input file.  Found the key '" + key
                                   + "' inside '" + c_state.name + "'.");
        }
      else if(parsing_B)
        {
          throw std::runtime_error("Invalid input file.  Found the key '" + key
                                   + "' inside '" + B_state.name + "'.");
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
          throw std::runtime_error("Invalid input file.  Found the key '" + key
                                   + "' inside the main object.");
        }
    }
  else
    {
      throw std::runtime_error("Found a key outside of the main object: "
                               + key);
    }
  return true;
}
