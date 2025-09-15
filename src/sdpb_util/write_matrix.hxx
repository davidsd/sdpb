#pragma once

#include <El.hpp>

template <class T>
void write_matrix(const El::Matrix<T> &matrix,
                  const std::filesystem::path &path)
{
  if(path.has_parent_path())
    std::filesystem::create_directories(path.parent_path());
  std::ofstream stream(path);
  El::Print(matrix,
            std::to_string(matrix.Height()) + " "
              + std::to_string(matrix.Width()),
            "\n", stream);
  stream << "\n";
  ASSERT(stream.good(), "Error when writing to: ", path);
}

// Code adapted from El::Print(const El::AbstractDistMatrix<T> &, ...)
// to work with with std::filesystem::path instead of ostream&
template <class T>
void write_distmatrix(const El::AbstractDistMatrix<T> &A,
                      const std::filesystem::path &path)
{
  if(A.ColStride() == 1 && A.RowStride() == 1)
    {
      if(A.CrossRank() == A.Root() && A.RedundantRank() == 0)
        {
          write_matrix(A.LockedMatrix(), path);
        }
    }
  else
    {
      El::DistMatrix<T, El::CIRC, El::CIRC> A_CIRC_CIRC(A);
      if(A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root())
        {
          write_matrix(A_CIRC_CIRC.LockedMatrix(), path);
        }
    }
}
