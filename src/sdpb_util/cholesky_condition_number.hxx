#pragma once

#include "assert.hxx"
#include <El.hpp>

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