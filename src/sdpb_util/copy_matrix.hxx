#pragma once

#include <El.hpp>

inline void copy_matrix(const El::Matrix<El::BigFloat> &source,
                        El::DistMatrix<El::BigFloat> &destination)
{
  destination.Resize(source.Height(), source.Width());
  for(int64_t row(0); row < destination.LocalHeight(); ++row)
    {
      int64_t global_row(destination.GlobalRow(row));
      for(int64_t column(0); column < destination.LocalWidth(); ++column)
        {
          int64_t global_column(destination.GlobalCol(column));
          destination.SetLocal(row, column, source(global_row, global_column));
        }
    }
}

inline void
copy_matrix(const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &source,
            El::Matrix<El::BigFloat> &destination)
{
  destination.Resize(source.LocalHeight(), source.LocalWidth());
  for(int64_t row(0); row < source.LocalHeight(); ++row)
    {
      for(int64_t column(0); column < source.LocalWidth(); ++column)
        {
          destination(row, column) = source.GetLocal(row, column);
        }
    }
}

inline void
copy_matrix(const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &source,
            El::DistMatrix<El::BigFloat> &destination)
{
  destination.Resize(source.LocalHeight(), source.LocalWidth());
  for(int64_t row(0); row < destination.LocalHeight(); ++row)
    {
      const int64_t global_row(destination.GlobalRow(row));
      for(int64_t column(0); column < destination.LocalWidth(); ++column)
        {
          const int64_t global_column(destination.GlobalCol(column));
          destination.SetLocal(row, column,
                               source.GetLocal(global_row, global_column));
        }
    }
}
