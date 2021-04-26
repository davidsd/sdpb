#pragma once

#include <El.hpp>

#include <boost/filesystem.hpp>

void read_text_block(El::DistMatrix<El::BigFloat> &block,
                     const boost::filesystem::path &block_path);
void read_text_block(El::DistMatrix<El::BigFloat> &block,
                     const boost::filesystem::path &checkpoint_directory,
                     const std::string &prefix, const size_t &block_index);
