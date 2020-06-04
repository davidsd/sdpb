#include "../Damped_Rational_State.hxx"

bool Damped_Rational_State::xml_on_characters(const xmlChar *characters,
                                              int length)
{
  if(inside)
    {
      constant_state.xml_on_characters(characters, length)
        || polynomial_state.xml_on_characters(characters, length)
        || base_state.xml_on_characters(characters, length);
    }
  return inside;
}
