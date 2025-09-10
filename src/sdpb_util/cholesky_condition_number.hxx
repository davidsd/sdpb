#pragma once

#include "assert.hxx"
#include <El.hpp>
#include <optional>

// Condition number for Cholesky matrix can be estimated as a squared ratio of maximal and minimal elements on its diagonal,
// see https://scicomp.stackexchange.com/questions/32762/cholmod-condition-number-estimate
template <class T>
T cholesky_condition_number(const El::AbstractDistMatrix<T> &matrix)
{
  ASSERT_EQUAL(matrix.Height(), matrix.Width());

  auto max_diagonal_element = El::limits::Min<T>();
  auto min_diagonal_element = El::limits::Max<T>();
  for(int i = 0; i < matrix.Height(); ++i)
    {
      if(matrix.IsLocal(i, i))
        {
          const int iLoc = matrix.LocalRow(i);
          const int jLoc = matrix.LocalCol(i);
          const auto &elt = matrix.GetLocalCRef(iLoc, jLoc);
          ASSERT(elt >= 0,
                 "Diagonal elements (i,i) of Cholesky matrix must be "
                 "non-negative!",
                 DEBUG_STRING(i), DEBUG_STRING(elt));
          max_diagonal_element = std::max(max_diagonal_element, elt);
          min_diagonal_element = std::min(min_diagonal_element, elt);
        }
    }
  max_diagonal_element = El::mpi::AllReduce(max_diagonal_element, El::mpi::MAX,
                                            matrix.DistComm());
  min_diagonal_element = El::mpi::AllReduce(min_diagonal_element, El::mpi::MIN,
                                            matrix.DistComm());

  auto ratio = max_diagonal_element / min_diagonal_element;
  return ratio * ratio;
}

// Returns max condition number among all blocks a distributed block-diagonal matrix
// and returns them on rank = 0.
// Returns {-1, -1} on other ranks
template <class TDistMatrix>
std::pair<int, El::BigFloat> max_block_cholesky_condition_number(
  const std::vector<TDistMatrix> &local_blocks, const size_t num_blocks,
  const std::function<size_t(size_t)> get_block_index)
{
  El::BigFloat zero(0);

  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> tttt = local_blocks.at(0);
  El::DistMatrix<El::BigFloat, El::CIRC, El::CIRC> wwww = local_blocks.at(0);

  std::vector<El::BigFloat> cond_numbers(num_blocks, zero);
  for(size_t i = 0; i < local_blocks.size(); ++i)
    {
      const auto block_index = get_block_index(i);
      ASSERT(block_index < num_blocks, DEBUG_STRING(block_index),
             DEBUG_STRING(num_blocks));
      cond_numbers.at(block_index)
        = cholesky_condition_number(local_blocks.at(i));
    }

  // Synchronize all at rank 0
  El::mpi::Reduce(cond_numbers.data(), num_blocks, El::mpi::MAX, 0,
                  El::mpi::COMM_WORLD);

  El::BigFloat max_cond_number = -1;
  int max_block_index = -1;
  if(El::mpi::Rank() == 0)
    {
      for(size_t i = 0; i < num_blocks; ++i)
        {
          if(cond_numbers.at(i) > max_cond_number)
            {
              max_cond_number = cond_numbers.at(i);
              max_block_index = i;
            }
        }
    }
  return {max_block_index, max_cond_number};
}
