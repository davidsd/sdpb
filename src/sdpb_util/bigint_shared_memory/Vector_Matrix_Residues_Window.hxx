#pragma once

#include "Array_Residues_Window.hxx"
#include "sdpb_util/assert.hxx"

// Residues for vector<Matrix>, arranged in a shared memory buffer as:
// matrices_residues[prime_index=0][matrix_index=0]
// matrices_residues[prime_index=0][matrix_index=1]
// ...
// matrices_residues[prime_index=0][matrix_index=num_matrices-1]
// matrices_residues[prime_index=1][matrix_index=0]
// ...
template <class T>
class Vector_Matrix_Residues_Window : public Array_Residues_Window<T>
{
public:
  const size_t num_matrices;
  // Indexed by (prime_index, matrix_index)
  std::vector<std::vector<El::Matrix<T>>> matrices_residues;

  Vector_Matrix_Residues_Window(Shared_Window_Array_View<T> window_view,
                                size_t num_primes,
                                const std::vector<size_t> &block_heights,
                                const std::vector<size_t> &block_widths)
      : Array_Residues_Window<T>(window_view, num_primes,
                                 num_elements(block_heights, block_widths)),
        num_matrices(block_heights.size()),
        matrices_residues(num_primes, std::vector<El::Matrix<T>>(num_matrices))
  {
    for(size_t prime_index = 0; prime_index < num_primes; ++prime_index)
      {
        size_t offset = 0;
        for(size_t i = 0; i < num_matrices; ++i)
          {
            const size_t height = block_heights.at(i);
            const size_t width = block_widths.at(i);
            T *data = this->data(prime_index) + offset;

            // Sanity checks
            ASSERT(offset <= this->prime_stride,
                   "Too many matrix elements, residues will overlap",
                   DEBUG_STRING(prime_index), DEBUG_STRING(i),
                   DEBUG_STRING(offset), DEBUG_STRING(this->prime_stride),
                   DEBUG_STRING(height), DEBUG_STRING(width));
            ASSERT(data + height * width <= this->data(0) + this->size(),
                   "Index out of Shared_Window_Array_View bounds");

            matrices_residues.at(prime_index)
              .at(i)
              .Attach(height, width, data, height);
            offset += height * width;
          }
      }
  }

  // Total number of matrix elements
  template <class TNumber>
  static TNumber num_elements(const std::vector<TNumber> &block_heights,
                              const std::vector<TNumber> &block_widths)
  {
    TNumber result(0);
    const auto num_blocks = block_heights.size();
    ASSERT_EQUAL(block_heights.size(), block_widths.size());
    for(size_t i = 0; i < num_blocks; ++i)
      {
        result += block_heights.at(i) * block_widths.at(i);
      }
    return result;
  }
};
