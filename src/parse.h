//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#ifndef SDPB_PARSE_H_
#define SDPB_PARSE_H_

//Tweak to allow Ubuntu-14.04/gcc-4.8.4 and similar environments to compile
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include "SDP.h"

using boost::filesystem::path;

SDP readBootstrapSDP(const vector<path> sdpFiles);

#endif  // SDPB_PARSE_H_
