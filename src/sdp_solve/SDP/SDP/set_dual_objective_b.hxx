#include <El.hpp>

void set_dual_objective_b(const std::vector<El::BigFloat> &objectives_input,
                          const El::Grid &grid,
                          El::DistMatrix<El::BigFloat> &dual_objective_b);
