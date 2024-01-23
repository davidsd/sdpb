#pragma once

#include "Polynomial_Vector_Matrix.hxx"

#include <El.hpp>

#include <vector>

// Result of JSON/Mathematica/xml input parsing.
// Format is described in SDPB Manual, eq. (3.1)
// NB: SDPB internals use format (2.2), which doesn't have normalization vector.
// Conversion from (3.1) to (2.2)
// takes place in pmp2sdp/Dual_Constraint_Group constructor,
// before writing to sdp.zip
struct Polynomial_Matrix_Program
{
  // vector a_i, i=0..N
  std::vector<El::BigFloat> objective;
  // normalization vector n_i, i=0..N
  std::vector<El::BigFloat> normalization;
  // Total number of matrices
  size_t num_matrices = 0;
  // In case of several processes,
  // each process owns only some matrices.
  std::vector<Polynomial_Vector_Matrix> matrices;
  // global index of matrices[i], lies in [0..num_matrices)
  std::vector<size_t> matrix_index_local_to_global;

  [[nodiscard]]
  Polynomial_Matrix_Program(std::vector<El::BigFloat> objective,
                            std::vector<El::BigFloat> normalization,
                            size_t num_matrices,
                            std::vector<Polynomial_Vector_Matrix> matrices,
                            std::vector<size_t> matrix_index_local_to_global)
      : objective(std::move(objective)),
        normalization(std::move(normalization)),
        num_matrices(num_matrices),
        matrices(std::move(matrices)),
        matrix_index_local_to_global(std::move(matrix_index_local_to_global))
  {
    // Validate
    if(this->num_matrices == 0)
      El::RuntimeError("Polynomial_Matrix_Program is empty!");
    if(this->objective.empty())
      El::RuntimeError("objective is empty!");
    if(this->objective.size() != this->normalization.size())
      {
        El::RuntimeError("objective.size=", this->objective.size(),
                         "should be equal to normalization.size=",
                         this->normalization.size());
      }
    if(this->matrices.size() > num_matrices)
      {
        El::RuntimeError("matrices.size=", this->matrices.size(),
                         " should not exceed", num_matrices);
      }
    if(this->matrices.size() != this->matrix_index_local_to_global.size())
      {
        El::RuntimeError(
          "matrices.size=", this->matrices.size(),
          " should be equal to matrix_index_local_to_global.size",
          this->matrix_index_local_to_global.size());
      }
    for(const size_t global_index : this->matrix_index_local_to_global)
      {
        if(global_index >= num_matrices)
          {
            El::RuntimeError(
              "Block index=", global_index,
              " should be less than num_matrices=", num_matrices);
          }
      }
    // TODO: we should also check that matrix indices from all ranks
    // are unique and cover [0, num_matrices) range.
    // This is checked indirectly in write_sdp().
  }

  // Allow move and prohibit copy

  Polynomial_Matrix_Program(const Polynomial_Matrix_Program &other) = delete;
  Polynomial_Matrix_Program(Polynomial_Matrix_Program &&other) noexcept
    = default;
  Polynomial_Matrix_Program &operator=(const Polynomial_Matrix_Program &other)
    = delete;
  Polynomial_Matrix_Program &
  operator=(Polynomial_Matrix_Program &&other) noexcept
    = default;
};
