#include "../Input_Parser.hxx"

void Input_Parser::on_characters(const Glib::ustring &characters)
{
  if(inside_sdp)
    {
      if(!objective_state.on_characters(characters))
        {
          polynomial_vector_matrices_state.on_characters(characters);
        }
    }
}
