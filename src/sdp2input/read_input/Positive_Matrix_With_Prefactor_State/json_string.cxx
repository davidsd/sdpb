#include "../Positive_Matrix_With_Prefactor_State.hxx"

void Positive_Matrix_With_Prefactor_State::json_string(const std::string &s)
{
  if(parsing_damped_rational)
    {
      damped_rational_state.json_string(s);
    }
  else if(parsing_polynomials)
    {
      polynomials_state.json_string(s);
    }
  else
    {
      throw std::runtime_error("Invalid input file.  Unexpected "
                               "string inside '"
                               + name + "': '" + s + "'");
    }
}
