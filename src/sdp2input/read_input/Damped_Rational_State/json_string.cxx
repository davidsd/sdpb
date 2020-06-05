#include "../Damped_Rational_State.hxx"

void Damped_Rational_State::json_string(const std::string &s)
{
  if(parsing_constant)
    {
      constant_state.json_string(s);
    }
  else if(parsing_base)
    {
      base_state.json_string(s);
    }
  else if(parsing_polynomial)
    {
      polynomial_state.json_string(s);
    }
  else
    {
      throw std::runtime_error("Invalid input file.  Unexpected string '" + s
                               + "' inside '" + name + "'");
    }
}
