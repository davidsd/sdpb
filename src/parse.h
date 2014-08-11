#ifndef SDP_BOOTSTRAP_PARSE_H_
#define SDP_BOOTSTRAP_PARSE_H_

#include "boost/filesystem.hpp"
#include "SDP.h"

using boost::filesystem::path;

SDP readBootstrapSDP(const path sdpFile);

#endif  // SDP_BOOTSTRAP_PARSE_H_
