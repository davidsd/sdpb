// Author(s): D. Simmons-Duffin, April-October 2013

#include <vector>
#include "parse.h"
#include "tinyxml2.h"

using tinyxml2::XMLDocument;
using tinyxml2::XMLElement;

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

Real parseReal(XMLElement *rElt) {
  return Real(rElt->GetText());
}

Polynomial parsePolynomial(XMLElement *polElt) {
  Polynomial p;
  p.coeffs = parseMany("coeff", parseReal, polElt);
  return p;
}

vector<Real> parseVector(XMLElement *vecElt) {
  return parseMany("coord", parseReal, vecElt);
}

DampedPartialFraction parseDampedPartialFraction(XMLElement *dprElt) {
  Real rho              = parseReal(dprElt->FirstChildElement("rho"));
  Polynomial pol        = parsePolynomial(dprElt->FirstChildElement("polynomial"));
  vector<Real> poles    = parseMany("pole", parseReal, dprElt->FirstChildElement("poles"));
  vector<Real> residues = parseMany("residue", parseReal, dprElt->FirstChildElement("residues"));
  return DampedPartialFraction(rho, pol, poles, residues);
}

DampedPartialFractionVector parseDampedPartialFractionVector(XMLElement *vElt) {
  DampedPartialFractionVector v;
  v.components = parseMany("dampedPartialFraction", parseDampedPartialFraction, vElt);
  v.xLeft  = 0;
  v.xRight = 40;
  return v;
}

SemiInfiniteProgram parseSemiInfiniteProgram(XMLElement *sipElt) {
  SemiInfiniteProgram sip;
  sip.norm   = parseVector(sipElt->FirstChildElement("norm")->FirstChildElement("vector"));
  sip.points = parseMany("vector", parseVector, sipElt->FirstChildElement("vectors"));
  sip.curves = parseMany("dampedPartialFractionVector",
                         parseDampedPartialFractionVector,
                         sipElt->FirstChildElement("curves"));
  sip.pointObj.resize(sip.points.size());
  return sip;
}

SemiInfiniteProgram readSemiInfiniteProgram(const char*sipFile) {
  XMLDocument doc;
  doc.LoadFile(sipFile);
  return parseSemiInfiniteProgram(doc.FirstChildElement("semiInfiniteProgram"));
}
