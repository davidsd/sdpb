#include "../Damped_Rational_State.hxx"

bool Damped_Rational_State::xml_on_end_element(const std::string &element_name)
{
  bool result(inside);
  if(inside)
    {
      if(element_name != "Symbol")
        {
          if(constant_state.xml_on_end_element(element_name))
            {
              finished_constant = !constant_state.inside;
              if(finished_constant)
                {
                  using namespace std;
                  swap(value.constant, constant_state.value);
                }
            }
          else if(polynomial_state.xml_on_end_element(element_name))
            {
              if(!polynomial_state.inside)
                {
                  using namespace std;
                  swap(value.poles, polynomial_state.value);
                }
            }
          else if(base_state.xml_on_end_element(element_name))
            {
              finished_base = !base_state.inside;
              if(finished_base)
                {
                  using namespace std;
                  swap(value.base, base_state.value);
                }
            }
          else
            {
              inside = (element_name != name);
            }
        }
    }
  return result;
}
