#pragma once

#include <El.hpp>

#include <vector>

void solve_LP(const El::Matrix<El::BigFloat> &A,
              const El::Matrix<El::BigFloat> &b,
              const El::Matrix<El::BigFloat> &c,
              std::vector<El::BigFloat> &weights);
