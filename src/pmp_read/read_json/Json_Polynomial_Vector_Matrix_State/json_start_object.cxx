#include "../Json_Polynomial_Vector_Matrix_State.hxx"

void Json_Polynomial_Vector_Matrix_State::json_start_object()
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
