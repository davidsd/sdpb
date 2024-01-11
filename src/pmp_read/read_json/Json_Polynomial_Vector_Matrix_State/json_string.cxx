#include "../Json_Polynomial_Vector_Matrix_State.hxx"

void Json_Polynomial_Vector_Matrix_State::json_string(const std::string &s)
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
