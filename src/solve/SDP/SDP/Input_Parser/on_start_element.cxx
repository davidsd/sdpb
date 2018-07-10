#include "../Input_Parser.hxx"

#include <stdexcept>

void Input_Parser::on_start_element(const Glib::ustring &element_name,
                                    const AttributeList &attributes)
{
  if(inside_sdp)
    {
      if(!objective_state.on_start_element(element_name, attributes)
         && !polynomial_vector_matrices_state.on_start_element(element_name,
                                                               attributes))
        {
          throw std::runtime_error(
            ("Invalid input file.  Expected '" + objective_state.name
             + "' or '" + polynomial_vector_matrices_state.name
             + "' inside 'sdp', but found '" + element_name + "'")
              .c_str());
        }
    }
  else if(glib_equals_string(element_name,sdp_name))
    {
      inside_sdp = true;
    }
  else
    {
      throw std::runtime_error(
        ("Invalid input file.  Expected 'sdp' but found '" + element_name
         + "'")
          .c_str());
    }
}
