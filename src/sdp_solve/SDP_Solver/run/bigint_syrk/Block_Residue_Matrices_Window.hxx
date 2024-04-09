#pragma once

#include "Residue_Matrices_Window.hxx"
#include "sdpb_util/assert.hxx"

// Same as Residue_Matrices_Window<T>,
// but each (tall) residue matrix (i.e. residues[prime_index])
// is split horizontally into blocks
// (i.e. block_residues[prime_index][0..num_blocks-1])
template <class T>
class Block_Residue_Matrices_Window : public Residue_Matrices_Window<T>
{
public:
  const size_t num_blocks;

  // block_residues[prime_index][block_index] = residue of block modulo prime
  // These matrices are views over residues[prime_index],
  // which is a tall matrix containing residues of each block,
  // stacked on top of each other:
  // block_residues[prime_index][0]
  // block_residues[prime_index][1]
  // ...
  // block_residues[prime_index][num_blocks-1]
  //
  // residues[prime_index] is a regular matrix attached to our memory window,
  // and block_residues[prime_index][block_index] is a view to its submatrix
  std::vector<std::vector<El::Matrix<T>>> block_residues;

  Block_Residue_Matrices_Window(El::mpi::Comm shared_memory_comm,
                                size_t num_primes, size_t num_blocks,
                                const std::vector<El::Int> &block_heights,
                                size_t block_width)
      : Residue_Matrices_Window<T>(shared_memory_comm, num_primes,
                                   Sum(block_heights), block_width),
        num_blocks(num_blocks),
        block_residues(num_primes, std::vector<El::Matrix<T>>(num_blocks))
  {
    for(size_t prime_index = 0; prime_index < num_primes; ++prime_index)
      {
        size_t block_start_row = 0;
        for(size_t block_index = 0; block_index < num_blocks; ++block_index)
          {
            size_t block_height = block_heights.at(block_index);
            El::Range<El::Int> I(block_start_row,
                                 block_start_row + block_height);
            El::Range<El::Int> J(0, this->width);
            El::View(block_residues.at(prime_index).at(block_index),
                     this->residues.at(prime_index), I, J);
            block_start_row += block_height;
          }
      }
    if(num_primes > 0 && num_blocks > 0)
      ASSERT_EQUAL(block_residues[0][0].Buffer(), this->residues[0].Buffer());
  }

private:
  static El::Int Sum(const std::vector<El::Int> &block_heights)
  {
    return std::accumulate(block_heights.begin(), block_heights.end(), 0);
  }
};
