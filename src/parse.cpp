//=======================================================================
// Copyright 2014 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#include <vector>
#include "types.h"
#include "parse.h"
#include "Polynomial.h"
#include "SDP.h"
#include "tinyxml2/tinyxml2.h"

using std::vector;
using tinyxml2::XMLDocument;
using tinyxml2::XMLElement;
using boost::filesystem::path;

template <class T>
vector<T> parseMany(const char *name, T(*parse)(XMLElement *),
                    XMLElement *elt) {
  XMLElement *e;
  vector<T> v;
  for (e = elt->FirstChildElement(name);
       e != NULL;
       e = e->NextSiblingElement(name)) {
    v.push_back(parse(e));
  }
  return v;
}

Real parseReal(XMLElement *xml) {
  return Real(xml->GetText());
}

int parseInt(XMLElement *xml) {
  return atoi(xml->GetText());
}

Vector parseVector(XMLElement *xml) {
  return parseMany("elt", parseReal, xml);
}

Matrix parseMatrix(XMLElement *xml) {
  Matrix m;
  m.rows     = parseInt(xml->FirstChildElement("rows"));
  m.cols     = parseInt(xml->FirstChildElement("cols"));
  m.elements = parseVector(xml->FirstChildElement("elements"));
  return m;
}

Polynomial parsePolynomial(XMLElement *xml) {
  Polynomial p;
  p.coefficients = parseMany("coeff", parseReal, xml);
  return p;
}

vector<Polynomial> parsePolynomialVector(XMLElement *xml) {
  return parseMany("polynomial", parsePolynomial, xml);
}

PolynomialVectorMatrix parsePolynomialVectorMatrix(XMLElement *xml) {
  PolynomialVectorMatrix m;
  m.rows           = parseInt(xml->FirstChildElement("rows"));
  m.cols           = parseInt(xml->FirstChildElement("cols"));
  m.elements       = parseMany("polynomialVector", parsePolynomialVector,
                               xml->FirstChildElement("elements"));
  m.samplePoints   = parseVector(xml->FirstChildElement("samplePoints"));
  m.sampleScalings = parseVector(xml->FirstChildElement("sampleScalings"));
  m.bilinearBasis  = parsePolynomialVector(xml->FirstChildElement("bilinearBasis"));
  return m;
}

SDP parseBootstrapSDP(XMLElement *xml) {
  return bootstrapSDP(parseVector(xml->FirstChildElement("objective")),
                      parseMany("polynomialVectorMatrix",
                                parsePolynomialVectorMatrix,
                                xml->FirstChildElement("polynomialVectorMatrices")));
}

SDP readBootstrapSDP(const path sdpFile) {
  XMLDocument doc;
  doc.LoadFile(sdpFile.c_str());
  return parseBootstrapSDP(doc.FirstChildElement("sdp"));
}

// Reading SDPs defined by DualConstraintGroups is not necessary --
// remove this functionality?

// DualConstraintGroup parseDualConstraintGroup(XMLElement *xml) {
//   DualConstraintGroup g;
//   g.dim                 = parseInt(xml->FirstChildElement("dim"));
//   g.degree              = parseInt(xml->FirstChildElement("degree"));
//   g.constraintMatrix    = parseMatrix(xml->FirstChildElement("constraintMatrix"));
//   g.constraintConstants = parseVector(xml->FirstChildElement("constraintConstants"));
//   g.bilinearBases       = parseMany("bilinearBasisMatrix",
//                                     parseMatrix,
//                                     xml->FirstChildElement("bilinearBases"));
//   return g;
// }

// SDP parseBootstrapSDP(XMLElement *xml) {
//   return sdpFromConstraintGroups(parseVector(xml->FirstChildElement("objective")),
//                                  parseReal(xml->FirstChildElement("objectiveConst")),
//                                  parseMany("dualConstraintGroup",
//                                            parseDualConstraintGroup,
//                                            xml->FirstChildElement("dualConstraintGroups")));
// }
