// See the manual for a description of the correct XML input format.

#include "read_xml.hxx"
#include "Input_Parser.hxx"

#include <filesystem>

namespace fs = std::filesystem;

namespace
{
  void start_element_callback(void *user_data, const xmlChar *name,
                              const xmlChar **)
  {
    Input_Parser *input_parser = static_cast<Input_Parser *>(user_data);
    input_parser->on_start_element(reinterpret_cast<const char *>(name));
  }

  void end_element_callback(void *user_data, const xmlChar *name)
  {
    Input_Parser *input_parser = static_cast<Input_Parser *>(user_data);
    input_parser->on_end_element(reinterpret_cast<const char *>(name));
  }

  void
  characters_callback(void *user_data, const xmlChar *characters, int length)
  {
    Input_Parser *input_parser = static_cast<Input_Parser *>(user_data);
    input_parser->on_characters(characters, length);
  }

  void warning_callback(void *, const char *msg, ...)
  {
    va_list args;
    va_start(args, msg);
    vprintf(msg, args);
    va_end(args);
  }

  void error_callback(void *, const char *msg, ...)
  {
    va_list args;
    va_start(args, msg);
    vprintf(msg, args);
    va_end(args);
    throw std::runtime_error("Invalid Input file");
  }
}

void read_xml(
  const std::filesystem::path &input_file,
  const std::function<bool(size_t matrix_index)> &should_parse_matrix,
  std::vector<El::BigFloat> &objective, size_t &num_matrices,
  std::map<size_t, Polynomial_Vector_Matrix> &polynomial_vector_matrices)
{
  LIBXML_TEST_VERSION;

  num_matrices = 0;
  auto process_matrix
    = [&should_parse_matrix, &num_matrices,
       &polynomial_vector_matrices](Polynomial_Vector_Matrix &&matrix) {
        // num_matrices equals to current matrix index
        auto index = num_matrices;
        if(should_parse_matrix(index))
          {
            polynomial_vector_matrices.emplace(index, matrix);
          }
        ++num_matrices;
      };
  Input_Parser input_parser(process_matrix);

  xmlSAXHandler xml_handlers;
  // This feels unclean.
  memset(&xml_handlers, 0, sizeof(xml_handlers));
  xml_handlers.startElement = start_element_callback;
  xml_handlers.endElement = end_element_callback;
  xml_handlers.characters = characters_callback;
  xml_handlers.warning = warning_callback;
  xml_handlers.error = error_callback;

  if(xmlSAXUserParseFile(&xml_handlers, &input_parser, input_file.c_str()) < 0)
    {
      throw std::runtime_error("Unable to parse input file: "
                               + input_file.string());
    }

  // Overwrite the objective with whatever is in the last file
  // that has an objective, but polynomial_vector_matrices get
  // appended.
  auto iterator(input_parser.objective_state.value.begin()),
    end(input_parser.objective_state.value.end());
  if(iterator != end)
    {
      objective.clear();
      objective.insert(objective.end(), iterator, end);
    }
}
