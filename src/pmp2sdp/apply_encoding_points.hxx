#pragma once

#include "Pmp2sdp_Parameters/Pmp2sdp_Parameters.hxx"
#include "pmp/Polynomial_Matrix_Program.hxx"

#include <El.hpp>

struct Encoding_Points_Report
{
  // Largest encoding point over the largest approximation point, over the
  // blocks on this rank
  double max_reach = 0;
};

// pmp2sdp --encodingPoints: set the encoding points and scalings of every
// block of pmp on this rank (no-op for --encodingPoints prefactorMatched).
Encoding_Points_Report
apply_encoding_points(Polynomial_Matrix_Program &pmp,
                      const Pmp2sdp_Parameters &parameters);

// Truncate to the given precision (the conversion runs at a higher one)
void truncate_precision(El::BigFloat &x, mp_bitcnt_t precision);
void truncate_precision(std::vector<El::BigFloat> &v, mp_bitcnt_t precision);
void truncate_precision(El::Matrix<El::BigFloat> &m, mp_bitcnt_t precision);
