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
  std::vector<El::BigFloat> objective;

  std::vector<Polynomial_Vector_Matrix> polynomialVectorMatrices;
  for(auto &sdp_file : sdp_files)
    {
      Input_Parser input_parser;
      input_parser.parse_file(sdp_file.string());

      std::swap(input_parser.objective_state.value, objective);
      std::swap(input_parser.polynomial_vector_matrices_state.value,
                polynomialVectorMatrices);
    }
  bootstrap(objective, polynomialVectorMatrices, *this);
}
