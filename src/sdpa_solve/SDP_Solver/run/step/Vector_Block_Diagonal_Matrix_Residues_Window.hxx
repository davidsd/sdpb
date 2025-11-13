#pragma once

#include "sdp_solve/Block_Matrix/Block_Diagonal_Matrix.hxx"
#include "sdpb_util/bigint_shared_memory/Vector_Matrix_Residues_Window.hxx"

namespace Sdpb::Sdpa
{
  // Residues for vector<Block_Diagonal_Matrix>
  template <class T>
  struct Vector_Block_Diagonal_Matrix_Residues_Window
      : Vector_Matrix_Residues_Window<T>
  {
    El::VerticalOrHorizontal stacking;

    Vector_Block_Diagonal_Matrix_Residues_Window(
      Shared_Window_Array_View<double> window_view, size_t num_primes,
      const std::vector<size_t> &block_dims, const size_t vector_size,
      const El::VerticalOrHorizontal stacking)
        : Vector_Matrix_Residues_Window<T>(
            window_view, num_primes,
            block_heights(block_dims, vector_size, stacking),
            block_widths(block_dims, vector_size, stacking)),
          stacking(stacking)
    {}

    El::Matrix<T>
    get_matrix_view(const size_t prime_index, const size_t block_index,
                    const size_t vec_index)
    {
      auto &stacked_matrix
        = this->matrices_residues.at(prime_index).at(block_index);
      const size_t dim = stacking == El::HORIZONTAL ? stacked_matrix.Height()
                                                    : stacked_matrix.Width();
      El::IR I(0, dim);
      El::IR J(dim * vec_index, dim * vec_index + dim);
      if(stacking == El::VERTICAL)
        {
          std::swap(I, J);
        }
      return stacked_matrix(I, J);
    }

  private:
    static std::vector<size_t>
    block_heights(std::vector<size_t> block_dims, const size_t vector_size,
                  const El::VerticalOrHorizontal stacking)
    {
      if(stacking == El::VERTICAL)
        {
          for(auto &dim : block_dims)
            dim *= vector_size;
        }
      return block_dims;
    }
    static std::vector<size_t>
    block_widths(std::vector<size_t> block_dims, const size_t vector_size,
                 const El::VerticalOrHorizontal stacking)
    {
      if(stacking == El::HORIZONTAL)
        {
          for(auto &dim : block_dims)
            dim *= vector_size;
        }
      return block_dims;
    }
  };
}