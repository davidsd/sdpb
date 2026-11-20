#pragma once

#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <filesystem>
#include <optional>
#include <vector>

// normalization.json of an SDP (a directory or a zip archive),
// at the precision of its decimals (at least SDPB's precision)
std::optional<std::vector<Boost_Float>>
read_normalization(const std::filesystem::path &sdp_path, Timers &timers);
