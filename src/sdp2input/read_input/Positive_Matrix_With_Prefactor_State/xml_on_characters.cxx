#include "../Positive_Matrix_With_Prefactor_State.hxx"

bool Positive_Matrix_With_Prefactor_State::xml_on_characters(
  const xmlChar *characters, int length)
{
  if(inside)
    {
      damped_rational_state.xml_on_characters(characters, length)
        || polynomials_state.xml_on_characters(characters, length);
    }
  return inside;
}
