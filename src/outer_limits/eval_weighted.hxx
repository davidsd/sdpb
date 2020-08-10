#pragma once

#include "../sdp_read.hxx"

El::BigFloat eval_weighted(const Positive_Matrix_With_Prefactor &matrix,
                           const El::BigFloat &x,
                           const std::vector<El::BigFloat> &weights);
