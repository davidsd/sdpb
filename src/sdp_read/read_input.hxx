#pragma once

#include "Positive_Matrix_With_Prefactor.hxx"

#include <boost/filesystem.hpp>

void read_input(const boost::filesystem::path &input_file,
                std::vector<El::BigFloat> &objectives,
                std::vector<El::BigFloat> &normalization,
                std::vector<Positive_Matrix_With_Prefactor> &matrices);

