#include "../Damped_Rational_State.hxx"

bool Damped_Rational_State::xml_on_start_element(
  const std::string &element_name)
{
  if(inside)
    {
      if(element_name != "Symbol"
         && (finished_constant
             || !constant_state.xml_on_start_element(element_name))
         && !polynomial_state.xml_on_start_element(element_name)
         && (finished_base || !base_state.xml_on_start_element(element_name)))
        {
          throw std::runtime_error(
            "Invalid input file.  Expected '" + constant_state.name + "' or '"
            + polynomial_state.name + "' inside DampedRational, but found '"
            + element_name + "'");
        }
    }
  else
    {
      inside = (element_name == name);
      if(inside)
        {
          finished_constant = false;
          finished_base = false;
        }
    }
  return inside;
}
