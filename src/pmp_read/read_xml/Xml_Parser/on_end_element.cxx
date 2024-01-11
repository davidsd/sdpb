#include "../Xml_Parser.hxx"

void Xml_Parser::on_end_element(const std::string &element_name)
{
  if(inside_sdp)
    {
      if(objective_state.xml_on_end_element(element_name))
        {
          finished_objective = !objective_state.inside;
        }
      else if(polynomial_vector_matrices_state.xml_on_end_element(
                element_name))
        {
          finished_polynomial_vector_matrices
            = !polynomial_vector_matrices_state.inside;
        }
      else if(element_name == sdp_name)
        {
          inside_sdp = false;
        }
    }
}
