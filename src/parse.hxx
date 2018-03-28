//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#pragma once

#include <boost/filesystem.hpp>
#include "SDP.hxx"

using boost::filesystem::path;

SDP readBootstrapSDP(const vector<path> sdpFiles);

