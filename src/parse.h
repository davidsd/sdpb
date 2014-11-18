//=======================================================================
// Copyright 2014 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#ifndef SDPB_PARSE_H_
#define SDPB_PARSE_H_

#include "boost/filesystem.hpp"
#include "SDP.h"

using boost::filesystem::path;

SDP readBootstrapSDP(const path sdpFile);

#endif  // SDPB_PARSE_H_
