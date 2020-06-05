#include "../Positive_Matrix_With_Prefactor_State.hxx"

void Positive_Matrix_With_Prefactor_State::json_start_object()
{
  if(inside)
    {
      if(parsing_damped_rational)
        {
          damped_rational_state.json_start_object();
        }
      else if(parsing_polynomials)
        {
          throw std::runtime_error(
            "Invalid input file.  Found an object inside '" + name + "."
            + polynomials_state.name + "'.");
        }
      else
        {
          throw std::runtime_error(
            "Invalid input file.  Unknown object inside '" + name + "'.");
        }
    }
  else
    {
      inside = true;
    }
}
