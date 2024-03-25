#pragma once

#include "assert.hxx"
#include "ostream/ostream_vector.hxx"

#include <El.hpp>
#include <iomanip>

// (TODO for debug only, don't merge) print diagonal elements and max/min ratio
inline void
print_cholesky_diagonal_info(const El::DistMatrix<El::BigFloat> &matrix,
                             const std::string &name)
{
  if(matrix.Height() == 0)
    return;

  ASSERT_EQUAL(matrix.Height(), matrix.Width());

  std::vector<double> diagonal(matrix.Height());
  for(int i = 0; i < matrix.Height(); ++i)
    {
      if(matrix.IsLocal(i, i))
        {
          const auto iLoc = matrix.LocalRow(i);
          const auto jLoc = matrix.LocalCol(i);
          diagonal.at(i) = static_cast<double>(matrix.GetLocal(iLoc, jLoc));
        }
    }

  const auto &comm = matrix.DistComm();
  El::mpi::Reduce(diagonal.data(), diagonal.size(), 0, comm);
  if(comm.Rank() == 0)
    {
      const double max_element
        = *std::max_element(diagonal.cbegin(), diagonal.cend());
      const double min_element
        = *std::min_element(diagonal.cbegin(), diagonal.cend());

      ASSERT(min_element > 0 && max_element >= min_element,
             DEBUG_STRING(min_element), DEBUG_STRING(max_element));

      std::ostringstream ss;
      ss << std::setprecision(3) << name
         << " diagonal: max/min ratio = " << max_element / min_element
         << ", elements = " << diagonal;
      El::Output(ss.str());
    }
}