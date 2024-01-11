#include "../Json_Polynomial_Vector_Matrix_State.hxx"

void Json_Polynomial_Vector_Matrix_State::json_key(const std::string &key)
{
  if(parsing_damped_rational)
    {
      damped_rational_state.json_key(key);
    }
  else if(parsing_polynomials)
    {
      throw std::runtime_error("Invalid input file.  Unexpected key '" + key
                               + "' inside '" + name + "."
                               + polynomials_state.name + "'.");
    }
  else if(key == damped_rational_state.name)
    {
      parsing_damped_rational = true;
    }
  else if(key == polynomials_state.name)
    {
      parsing_polynomials = true;
    }
  else
    {
      throw std::runtime_error("Invalid input file.  Unexpected key '" + key
                               + "' inside '" + name + "'");
    }
}
