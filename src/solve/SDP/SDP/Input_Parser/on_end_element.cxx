#include "../Input_Parser.hxx"

void Input_Parser::on_end_element(const Glib::ustring &element_name)
{
  if(inside_sdp)
    {
      if(objective_state.on_end_element(element_name))
        {
          finished_objective=true;
        }
      else if(polynomial_vector_matrices_state.on_end_element(element_name))
        {
          finished_polynomial_vector_matrices=true;
        }
      else if(glib_equals_string(element_name,sdp_name))
        {
          inside_sdp = false;
        }
    }
}
