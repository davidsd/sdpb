#include <vector>
#include "types.h"
#include "parse.h"
#include "tinyxml2.h"
#include "SDP.h"

using std::vector;
using tinyxml2::XMLDocument;
using tinyxml2::XMLElement;
using boost::filesystem::path;

template <class T>
vector<T> parseMany(const char *name, T(*parse)(XMLElement *), XMLElement *elt) {
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

SampledMatrixPolynomial parseSampledMatrixPolynomial(XMLElement *xml) {
  SampledMatrixPolynomial s;
  s.dim                 = parseInt(xml->FirstChildElement("dim"));
  s.degree              = parseInt(xml->FirstChildElement("degree"));
  s.constraintMatrix    = parseMatrix(xml->FirstChildElement("constraintMatrix"));
  s.constraintConstants = parseVector(xml->FirstChildElement("constraintConstants"));
  s.bilinearBases       = parseMany("bilinearBasisMatrix", parseMatrix, xml->FirstChildElement("bilinearBases"));
  return s;
}


SDP parseBootstrapSDP(XMLElement *xml) {
  return bootstrapSDP(parseVector(xml->FirstChildElement("objective")),
                      parseReal(xml->FirstChildElement("objectiveConst")),
                      parseMany("sampledMatrixPolynomial",
                                parseSampledMatrixPolynomial,
                                xml->FirstChildElement("sampledPositiveMatrices")));
}

SDP readBootstrapSDP(const path sdpFile) {
  XMLDocument doc;
  doc.LoadFile(sdpFile.c_str());
  return parseBootstrapSDP(doc.FirstChildElement("sdp"));
}
