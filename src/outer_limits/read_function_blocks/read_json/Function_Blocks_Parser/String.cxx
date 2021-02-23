#include "../Function_Blocks_Parser.hxx"

bool Function_Blocks_Parser::String(const Ch *str, rapidjson::SizeType length,
                                    bool)
{
  std::string s(str, length);
  if(inside)
    {
      if(parsing_objective)
        {
          objective_state.json_string(s);
        }
      else if(parsing_normalization)
        {
          normalization_state.json_string(s);
        }
      else if(parsing_functions)
        {
          functions_state.json_string(s);
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
      throw std::runtime_error("Found a string outside of the SDP: '" + s
                               + "'");
    }
  return true;
}
