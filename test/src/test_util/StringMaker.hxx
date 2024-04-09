#pragma once

#include "catch2/catch_amalgamated.hpp"
#include <El.hpp>

// Custom functions to log

namespace Test_Util
{
  template <class TMatrix>
  std::string to_string(const TMatrix &matrix, const std::string &label)
  {
    std::ostringstream os;
    os << El::DimsString(matrix, label);
    os << "\n[\n";
    for(int i = 0; i < matrix.Height(); ++i)
      {
        for(int j = 0; j < matrix.Width(); ++j)
          {
            os << " " << matrix.Get(i, j);
          }
        os << "\n";
      }
    os << "]";
    return os.str();
  }
}

// Catch::StringMaker is used to print variables in CAPTURE() macro.
namespace Catch
{
  template <> struct StringMaker<El::Matrix<double>>
  {
    static std::string convert(const El::Matrix<double> &matrix)
    {
      return Test_Util::to_string(matrix, "Matrix<double>");
    }
  };

  template <> struct StringMaker<El::Matrix<El::BigFloat>>
  {
    static std::string convert(const El::Matrix<El::BigFloat> &matrix)
    {
      return Test_Util::to_string(matrix, "Matrix<BigFloat>");
    }
  };

  template <> struct StringMaker<El::DistMatrix<El::BigFloat>>
  {
    static std::string convert(const El::DistMatrix<El::BigFloat> &matrix)
    {
      return Test_Util::to_string(matrix, "DistMatrix<BigFloat>");
    }
  };
}
