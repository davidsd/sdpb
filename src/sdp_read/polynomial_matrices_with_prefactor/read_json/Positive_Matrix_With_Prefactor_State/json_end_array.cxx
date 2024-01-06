#include "../Positive_Matrix_With_Prefactor_State.hxx"

void Positive_Matrix_With_Prefactor_State::json_end_array()
{
  if(parsing_damped_rational)
    {
      damped_rational_state.json_end_array();
    }
  else if(parsing_polynomials)
    {
      polynomials_state.json_end_array();
      parsing_polynomials=polynomials_state.inside;
      if(!parsing_polynomials)
        {
          swap(value.polynomials,polynomials_state.value);
        }
    }
  else
    {
      throw std::runtime_error(
        "Invalid input file.  Unknown array end inside '" + name + "'");
    }
}
