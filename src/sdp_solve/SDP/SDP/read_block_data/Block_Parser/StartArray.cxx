#include "Block_Parser.hxx"

bool Block_Parser::StartArray()
{
  if(inside)
    {
      if(parsing_c)
        {
          c_state.json_start_array();
        }
      else if(parsing_B)
        {
          B_state.json_start_array();
        }
      else if(parsing_bilinear_bases_even)
        {
          bilinear_bases_even_state.json_start_array();
        }
      else if(parsing_bilinear_bases_odd)
        {
          bilinear_bases_odd_state.json_start_array();
        }
      else
        {
          throw std::runtime_error(
            "Invalid input file.  Found an array not inside " + c_state.name
            + ", " + B_state.name + ", " + bilinear_bases_even_state.name
            + ", " + bilinear_bases_odd_state.name + ".");
        }
    }
  else
    {
      throw std::runtime_error("Found an array outside of the main object.");
    }
  return true;
}
