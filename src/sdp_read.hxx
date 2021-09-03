#pragma once

#include "sdp_read/read_input.hxx"
#include "sdp_read/read_pvm_input.hxx"

std::vector<boost::filesystem::path>
read_file_list(const boost::filesystem::path &input_file);
