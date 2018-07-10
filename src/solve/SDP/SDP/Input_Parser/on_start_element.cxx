#include "../Input_Parser.hxx"

#include <stdexcept>

void Input_Parser::on_start_element(const Glib::ustring &name,
                                    const AttributeList &attributes)
{
  if(name == "sdp")
    {
      inside_sdp = true;
    }
  else if(inside_sdp)
    {
      if(!objective_state.on_start_element(name, attributes)
         && !polynomial_vector_matrices_state.on_start_element(name,
                                                               attributes))
        {
          throw std::runtime_error(("Invalid input file.  Expected '"
          + objective_state.name + "' or '"
          + polynomial_vector_matrices_state.name
                              + "' inside 'sdp', but found '" + name + "'").c_str());
        }
    }
  else
    {
      throw std::runtime_error(
        ("Invalid input file.  Expected 'sdp' but found '" + name + "'").c_str());
    }
}
