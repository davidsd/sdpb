//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#ifndef SDPB_SERIALIZE_H_
#define SDPB_SERIALIZE_H_

#include <string>
#include <vector>
#include <sstream>
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/base_object.hpp"
#include "boost/serialization/utility.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/string.hpp"
#include "boost/serialization/split_free.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/archive/text_oarchive.hpp"
//Tweak to allow Ubuntu-14.04/gcc-4.8.4 and similar environments to compile
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include "boost/filesystem/fstream.hpp"
#include "types.h"
#include "Vector.h"
#include "Matrix.h"
#include "BlockDiagonalMatrix.h"

using boost::filesystem::path;
using boost::archive::text_iarchive;
using std::vector;

namespace boost {
  namespace serialization {

    template<class Archive>
    void load(Archive& ar, Real& f, unsigned int version) {
      std::string s;
      ar & s;
      f = Real(s.c_str());
    }

    template<class Archive>
    void save(Archive& ar, Real const& f, unsigned int version) {
      std::ostringstream os;
      os.precision(f.get_prec());
      os << f;
      std::string s = os.str();
      ar & s;
    }

    template<class Archive>
    void serialize(Archive& ar, Matrix& m, const unsigned int version) {
      ar & m.rows;
      ar & m.cols;
      ar & m.elements;
    }

    template<class Archive>
    void serialize(Archive& ar, BlockDiagonalMatrix& m, const unsigned int version) {
      ar & m.dim;
      ar & m.blocks;
      ar & m.blockStartIndices;
    }

    template<class Archive>
    void serializeSDPSolverState(Archive& ar, Vector &x, BlockDiagonalMatrix &X, Vector &y, BlockDiagonalMatrix &Y) {
      ar & x;
      ar & X;
      ar & y;
      ar & Y;
    }

  }  // namespace serialization
}  // namespace boost


BOOST_SERIALIZATION_SPLIT_FREE(Real)
BOOST_CLASS_VERSION(Real, 0)
BOOST_CLASS_TRACKING(Matrix,              boost::serialization::track_never)
BOOST_CLASS_TRACKING(BlockDiagonalMatrix, boost::serialization::track_never)

#endif  // SDPB_SERIALIZE_H_
