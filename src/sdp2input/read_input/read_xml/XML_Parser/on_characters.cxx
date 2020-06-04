#include "../XML_Parser.hxx"

void XML_Parser::on_characters(const xmlChar *characters, int length)
{
  if(inside_expression && inside_sdp)
    {
      objective_state.xml_on_characters(characters, length)
        || normalization_state.xml_on_characters(characters, length)
        || positive_matrices_with_prefactor_state.xml_on_characters(characters,
                                                                    length);
    }
}
