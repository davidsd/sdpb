//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

// See the manual for a description of the correct XML input format.

#include "Input_Parser.hxx"
#include "../../SDP.hxx"

#include <boost/filesystem.hpp>
#include <boost/property_tree/xml_parser.hpp>

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
}

void bootstrap(const std::vector<El::BigFloat> &objective,
               const std::vector<Polynomial_Vector_Matrix> &polVectorMatrices,
               SDP &sdp);

Polynomial_Vector_Matrix
parse_polynomial_vector_matrix(const boost::property_tree::ptree &tree);

SDP::SDP(const std::vector<boost::filesystem::path> &sdp_files)
    // FIXME: This assigns one core per block.  We may want to do
    // something more sophisticated for larger blocks.
    : grid(El::mpi::COMM_SELF)
{
  LIBXML_TEST_VERSION;

  std::vector<El::BigFloat> objective;

  std::vector<Polynomial_Vector_Matrix> polynomialVectorMatrices;
  for(auto &sdp_file : sdp_files)
    {
      Input_Parser input_parser;

      xmlSAXHandler xml_handlers = {0};
      xml_handlers.startElement = start_element_callback;
      xml_handlers.endElement = end_element_callback;
      xml_handlers.characters = characters_callback;

      if(xmlSAXUserParseFile(&xml_handlers, &input_parser, sdp_file.c_str())
         < 0)
        {
          throw std::runtime_error("Ill-formed input file: "
                                   + sdp_file.string());
        }

      std::swap(input_parser.objective_state.value, objective);
      std::swap(input_parser.polynomial_vector_matrices_state.value,
                polynomialVectorMatrices);
    }
  bootstrap(objective, polynomialVectorMatrices, *this);
}
