//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


// Currently, SDPB uses tinyxml2 to parse an input file into an XML
// tree, which is stored entirely in memory.  This tree is then
// transformed into the appropriate data structures using the parse
// functions below.  In the future, this could be made more memory
// efficient by avoiding building the XML tree in memory.

// See the manual for a description of the correct XML input format.

#include "parse_append_many.hxx"

#include "../Polynomial.hxx"
#include "../SDP.hxx"

#include <boost/filesystem.hpp>
#include <boost/property_tree/xml_parser.hpp>

std::vector<Real>
parse_vector(const boost::property_tree::ptree &tree);

Polynomial_Vector_Matrix parse_polynomial_vector_matrix(const boost::property_tree::ptree &tree);

SDP read_bootstrap_sdp(const std::vector<boost::filesystem::path> sdp_files)
{
  Vector objective;
  std::vector<Polynomial_Vector_Matrix> polynomialVectorMatrices;
  for (auto &sdp_file: sdp_files)
    {
      {
      boost::property_tree::ptree tree;
      boost::property_tree::read_xml (sdp_file.string(), tree);
      
      const auto sdp = tree.get_child ("sdp");
      auto objective_iterator (sdp.find("objective"));
      /// boost::property_tree uses not_found() instead of end() :(
      if(objective_iterator!=sdp.not_found())
        { objective=parse_vector(objective_iterator->second); }

      auto polynomialVectorMatrices_iterator (sdp.find("polynomialVectorMatrices"));
      if(polynomialVectorMatrices_iterator!=sdp.not_found())
        {
          std::function<Polynomial_Vector_Matrix(const boost::property_tree::ptree &)>
            p (parse_polynomial_vector_matrix);
          parse_append_many("polynomialVectorMatrix",
                            p,
                            polynomialVectorMatrices_iterator->second,
                            polynomialVectorMatrices);
                            
        }
      }
    }
  return bootstrapSDP(objective,polynomialVectorMatrices);
}
