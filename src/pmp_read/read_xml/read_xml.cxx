// See the manual for a description of the correct XML input format.

#include "Xml_Parser.hxx"
#include "pmp_read/PMP_File_Parse_Result.hxx"

#include <filesystem>

namespace fs = std::filesystem;

namespace
{
  void start_element_callback(void *user_data, const xmlChar *name,
                              const xmlChar **)
  {
    Xml_Parser *input_parser = static_cast<Xml_Parser *>(user_data);
    input_parser->on_start_element(reinterpret_cast<const char *>(name));
  }

  void end_element_callback(void *user_data, const xmlChar *name)
  {
    Xml_Parser *input_parser = static_cast<Xml_Parser *>(user_data);
    input_parser->on_end_element(reinterpret_cast<const char *>(name));
  }

  void
  characters_callback(void *user_data, const xmlChar *characters, int length)
  {
    Xml_Parser *input_parser = static_cast<Xml_Parser *>(user_data);
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

PMP_File_Parse_Result
read_xml(const std::filesystem::path &input_file,
         const std::function<bool(size_t matrix_index)> &should_parse_matrix)
{
  LIBXML_TEST_VERSION;

  PMP_File_Parse_Result result;

  result.num_matrices = 0;
  // TODO this is ugly, we should keep track of matrix index inside parser
  // and skip matrices that we don't need.
  // NB: should be called only once for each matrix index,
  // otherwise num_matrices will be incorrect!
  // Currently it is called in Xml_Polynomial_Vector_Matrix_State.xml_on_end_element()
  auto process_matrix
    = [&should_parse_matrix, &result](
        Xml_Polynomial_Vector_Matrix_State &matrix_state) {
        // num_matrices equals to current matrix index
        auto index = result.num_matrices;
        if(should_parse_matrix(index))
          {
            El::Matrix<std::vector<Polynomial>> poly_vectors(
              matrix_state.rows, matrix_state.cols);

            auto elt = matrix_state.elements_state.value.begin();
            for(int i = 0; i < poly_vectors.Height(); ++i)
              for(int j = 0; j < poly_vectors.Width(); ++j)
                {
                  swap(poly_vectors(i, j), *elt++);
                }

            std::optional<Damped_Rational> prefactor = std::nullopt;
            auto &sample_points = matrix_state.sample_points_state.value;
            auto &sample_scalings = matrix_state.sample_scalings_state.value;
            std::vector<Polynomial> bilinear_basis;
            swap(bilinear_basis, matrix_state.bilinear_basis_state.value);

            result.parsed_matrices.emplace(
              index,
              Polynomial_Vector_Matrix(
                std::move(poly_vectors), prefactor, std::move(sample_points),
                std::move(sample_scalings), std::move(bilinear_basis)));
          }
        ++result.num_matrices;
      };
  Xml_Parser input_parser(process_matrix);

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
      result.objective.clear();
      result.objective.insert(result.objective.end(), iterator, end);
    }
  return result;
}
