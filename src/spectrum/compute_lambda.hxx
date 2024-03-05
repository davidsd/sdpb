#pragma once

#include "Zero.hxx"
#include "pmp_read/pmp_read.hxx"

void compute_lambda(const std::vector<El::BigFloat> &sample_points,
                    const std::vector<El::BigFloat> &sample_scalings,
                    const size_t &num_rows, const El::Matrix<El::BigFloat> &x,
                    const std::vector<El::BigFloat> &zero_vector,
                    std::vector<Zero> &zeros, El::BigFloat &error);
