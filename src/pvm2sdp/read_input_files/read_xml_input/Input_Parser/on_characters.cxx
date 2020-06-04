#include "../Input_Parser.hxx"

void Input_Parser::on_characters(const xmlChar *characters, int length)
{
  if(inside_sdp)
    {
      if(!objective_state.xml_on_characters(characters, length))
        {
          polynomial_vector_matrices_state.xml_on_characters(characters,
                                                             length);
        }
    }
}
