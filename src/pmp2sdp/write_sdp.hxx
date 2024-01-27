#pragma once

#include "Block_File_Format.hxx"
#include "Output_SDP/Output_SDP.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <filesystem>

void write_sdp(const std::filesystem::path &output_path, const Output_SDP &sdp,
               Block_File_Format block_file_format, bool zip, Timers &timers,
               bool debug);
