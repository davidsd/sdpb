#include "XML_Parser.hxx"
#include "../../sdp_convert.hxx"

#include <boost/filesystem.hpp>

void read_mathematica(const boost::filesystem::path &input_path,
                      std::vector<El::BigFloat> &objectives,
                      std::vector<El::BigFloat> &normalization,
                      std::vector<Positive_Matrix_With_Prefactor> &matrices);

namespace
{
  void start_element_callback(void *user_data, const xmlChar *name,
                              const xmlChar **)
  {
    XML_Parser *input_parser = static_cast<XML_Parser *>(user_data);
    input_parser->on_start_element(reinterpret_cast<const char *>(name));
  }

  void end_element_callback(void *user_data, const xmlChar *name)
  {
    XML_Parser *input_parser = static_cast<XML_Parser *>(user_data);
    input_parser->on_end_element(reinterpret_cast<const char *>(name));
  }

  void
  characters_callback(void *user_data, const xmlChar *characters, int length)
  {
    XML_Parser *input_parser = static_cast<XML_Parser *>(user_data);
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

void read_input(const boost::filesystem::path &input_file,
                std::vector<El::BigFloat> &objectives,
                std::vector<El::BigFloat> &normalization,
                std::vector<Positive_Matrix_With_Prefactor> &matrices)
{
  if(input_file.extension() == ".nsv")
    {
      for(auto &filename : read_file_list(input_file))
        {
          if(!filename.empty())
            {
              read_input(filename, objectives, normalization, matrices);
            }
        }
    }
  else if(input_file.extension() == ".xml")
    {
      LIBXML_TEST_VERSION;

      xmlSAXHandler xml_handlers;
      // This feels unclean.
      memset(&xml_handlers, 0, sizeof(xml_handlers));
      xml_handlers.startElement = start_element_callback;
      xml_handlers.endElement = end_element_callback;
      xml_handlers.characters = characters_callback;
      xml_handlers.warning = warning_callback;
      xml_handlers.error = error_callback;

      XML_Parser input_parser;
      if(xmlSAXUserParseFile(&xml_handlers, &input_parser, input_file.c_str())
         < 0)
        {
          throw std::runtime_error("Unable to parse input file: "
                                   + input_file.string());
        }

      if(!input_parser.objective_state.value.empty())
        {
          std::swap(objectives, input_parser.objective_state.value);
        }
      if(!input_parser.normalization_state.value.empty())
        {
          std::swap(normalization, input_parser.normalization_state.value);
        }
      size_t offset(matrices.size());
      auto &temp_matrices(
        input_parser.positive_matrices_with_prefactor_state.value);
      matrices.resize(matrices.size() + temp_matrices.size());
      for(size_t index = 0; index < temp_matrices.size(); ++index)
        {
          std::swap(matrices[offset + index], temp_matrices[index]);
        }
    }
  else
    {
      read_mathematica(input_file, objectives, normalization, matrices);
    }
}
