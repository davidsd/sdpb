#include "../Positive_Matrix_With_Prefactor_State.hxx"

void Positive_Matrix_With_Prefactor_State::json_key(const std::string &key)
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
