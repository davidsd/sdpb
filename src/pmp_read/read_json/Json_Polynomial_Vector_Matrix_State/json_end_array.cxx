#include "../Json_Polynomial_Vector_Matrix_State.hxx"

void Json_Polynomial_Vector_Matrix_State::json_end_array()
{
  if(parsing_damped_rational)
    {
      damped_rational_state.json_end_array();
    }
  else if(parsing_polynomials)
    {
      polynomials_state.json_end_array();
      parsing_polynomials=polynomials_state.inside;
    }
  else
    {
      throw std::runtime_error(
        "Invalid input file.  Unknown array end inside '" + name + "'");
    }
}
