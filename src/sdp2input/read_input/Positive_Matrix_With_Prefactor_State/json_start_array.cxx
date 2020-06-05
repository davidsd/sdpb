#include "../Positive_Matrix_With_Prefactor_State.hxx"

void Positive_Matrix_With_Prefactor_State::json_start_array()
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
