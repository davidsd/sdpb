#include "../Input_Parser.hxx"

void Input_Parser::on_end_element(const Glib::ustring &name)
{
  if(inside_sdp)
    {
      if(name == "sdp")
        {
          inside_sdp = false;
        }
      else if(objective_state.on_end_element(name))
        {
          finished_objective=true;
        }
      else if(polynomial_vector_matrices_state.on_end_element(name))
        {
          finished_polynomial_vector_matrices=true;
        }
    }
}
