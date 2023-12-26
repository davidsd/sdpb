#include "Block_Parser.hxx"

bool Block_Parser::String(const Ch *str, rapidjson::SizeType length, bool)
{
  std::string s(str, length);
  if(inside)
    {
      if(parsing_c)
        {
          c_state.json_string(s);
        }
      else if(parsing_B)
        {
          B_state.json_string(s);
        }
      else if(parsing_bilinear_bases_even)
        {
          bilinear_bases_even_state.json_string(s);
        }
      else if(parsing_bilinear_bases_odd)
        {
          bilinear_bases_odd_state.json_string(s);
        }
      else
        {
          throw std::runtime_error(
            "Invalid input file.  Unexpected string in the main object: '" + s
            + "'");
        }
    }
  else
    {
      throw std::runtime_error("Found a string outside of the main object: '"
                               + s + "'");
    }
  return true;
}
