#include "../Positive_Matrix_With_Prefactor_State.hxx"

void Positive_Matrix_With_Prefactor_State::json_end_object()
{
  if(parsing_damped_rational)
    {
      damped_rational_state.json_end_object();
      parsing_damped_rational=damped_rational_state.inside;
      if(!parsing_damped_rational)
        {
          swap(value.damped_rational,damped_rational_state.value);
        }
    }
  else if(parsing_polynomials)
    {
      throw std::runtime_error(
        "Invalid input file.  Found an object end inside '" + name + "."
        + polynomials_state.name + "'");
    }
  else
    {
      inside = false;
    }
}
