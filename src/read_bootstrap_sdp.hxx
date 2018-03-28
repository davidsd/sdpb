//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#pragma once

#include "SDP.hxx"
#include <boost/filesystem.hpp>

SDP read_bootstrap_sdp(const std::vector<boost::filesystem::path> sdpFiles);

