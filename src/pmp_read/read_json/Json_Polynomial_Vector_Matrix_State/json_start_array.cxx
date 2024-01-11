#include "../Json_Polynomial_Vector_Matrix_State.hxx"

void Json_Polynomial_Vector_Matrix_State::json_start_array()
{
  if(parsing_damped_rational)
    {
      damped_rational_state.json_start_array();
    }
  else if(parsing_polynomials)
    {
      polynomials_state.json_start_array();
    }
  else
    {
      throw std::runtime_error("Invalid input file.  Unknown array inside "
                               "'PositiveMatrixWithPrefactorArray'");
    }
}
