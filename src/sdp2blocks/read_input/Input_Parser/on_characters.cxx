#include "../Input_Parser.hxx"

void Input_Parser::on_characters(const xmlChar *characters, int length)
{
  if(inside_expression)
    {
      if(inside_sdp)
        {
          objective_state.on_characters(characters, length);
        }
    }
}
