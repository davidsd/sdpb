#include "../Xml_Parser.hxx"

#include <stdexcept>

void Xml_Parser::on_start_element(const std::string &element_name)
{
  if(inside_sdp)
    {
      if(!objective_state.xml_on_start_element(element_name)
         && !polynomial_vector_matrices_state.xml_on_start_element(
           element_name))
        {
          RUNTIME_ERROR("Invalid input file. Expected '", objective_state.name,
                        "' or '", polynomial_vector_matrices_state.name,
                        "' inside 'sdp', but found '", element_name, "'");
        }
    }
  else if(element_name == sdp_name)
    {
      inside_sdp = true;
    }
  else
    {
      RUNTIME_ERROR("Invalid input file. Expected 'sdp' but found '"
                    + element_name + "'");
    }
}
