#include "../Input_Parser.hxx"

void Input_Parser::on_characters(const xmlChar *characters, int length)
{
  // std::cout << "inside: " << std::string(reinterpret_cast<const char*>(characters),length) << "\n";
  
  if(inside_expression && inside_sdp)
    {
      objective_state.on_characters(characters, length)
        || normalization_state.on_characters(characters, length)
        || positive_matrix_with_prefactor_state.on_characters(characters, length);
    }
}
