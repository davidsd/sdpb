#include "../Positive_Matrix_With_Prefactor_State.hxx"

bool Positive_Matrix_With_Prefactor_State::xml_on_start_element(
  const std::string &element_name)
{
  if(inside)
    {
      if(element_name != "Symbol"
         && (finished_damped_rational
             || !damped_rational_state.xml_on_start_element(element_name))
         && !polynomials_state.xml_on_start_element(element_name))
        {
          throw std::runtime_error(
            "Invalid input file.  Expected 'Function' inside "
            "PositiveMatrixWithPrefactor, but found '"
            + element_name + "'");
        }
    }
  else
    {
      inside = (element_name == name);
      if(inside)
        {
          finished_damped_rational = false;
        }
    }
  return inside;
}
