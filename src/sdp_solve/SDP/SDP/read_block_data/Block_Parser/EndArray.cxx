#include "Block_Parser.hxx"

bool Block_Parser::EndArray(rapidjson::SizeType)
{
  if(inside)
    {
      if(parsing_c)
        {
          c_state.json_end_array();
          parsing_c = c_state.inside;
        }
      else if(parsing_B)
        {
          B_state.json_end_array();
          parsing_B = B_state.inside;
        }
      else if(parsing_bilinear_bases_even)
        {
          bilinear_bases_even_state.json_end_array();
          parsing_bilinear_bases_even = bilinear_bases_even_state.inside;
        }
      else if(parsing_bilinear_bases_odd)
        {
          bilinear_bases_odd_state.json_end_array();
          parsing_bilinear_bases_odd = bilinear_bases_odd_state.inside;
        }
      else
        {
          throw std::runtime_error(
            "INTERNAL ERROR: Ending an array inside the main object.");
        }
    }
  else
    {
      throw std::runtime_error(
        "Invalid input file.  Ending an array outside of the main object.");
    }
  return true;
}
