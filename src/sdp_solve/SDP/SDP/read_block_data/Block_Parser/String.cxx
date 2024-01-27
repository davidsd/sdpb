#include "Block_Parser.hxx"
#include "sdpb_util/assert.hxx"

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
          RUNTIME_ERROR(
            "Invalid input file. Unexpected string in the main object: '", s,
            "'");
        }
    }
  else
    {
      RUNTIME_ERROR("Found a string outside of the main object: '", s, "'");
    }
  return true;
}
