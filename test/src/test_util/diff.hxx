#pragma once

#include "pmp/Damped_Rational.hxx"
#include "pmp/PVM_Info.hxx"

#include <catch2/catch_amalgamated.hpp>
#include <El.hpp>
#include <boost/noncopyable.hpp>
#include <filesystem>

#define DIFF(a, b)                                                            \
  {                                                                           \
    INFO("DIFF(" << #a << ", " << #b << ")");                                 \
    diff(a, b);                                                               \
  }

// DIFF up to a given binary precision 'prec'
#define DIFF_PREC(a, b, prec)                                                 \
  {                                                                           \
    Test_Util::REQUIRE_Equal::Diff_Precision p(prec);                         \
    INFO("DIFF(" << #a << ", " << #b << ")");                                 \
    diff(a, b);                                                               \
  }

// Functions that REQUIRE equality of given files or folders.
// NB: we call REQUIRE() inside these functions,
// because otherwise useful CAPTURE() information will be lost
namespace Test_Util::REQUIRE_Equal
{
  inline int diff_precision;

  // RAII wrapper allowing to set diff_precision temporarily
  struct Diff_Precision : boost::noncopyable
  {
    explicit Diff_Precision(int precision)
    {
      old_precision = diff_precision;
      diff_precision = precision;
    }
    virtual ~Diff_Precision() { diff_precision = old_precision; }

  private:
    int old_precision;
  };

  template <class T> void diff(const T &a, const T &b)
  {
    REQUIRE(a == b);
  }

  template <> inline void diff(const El::BigFloat &a, const El::BigFloat &b)
  {
    if(a == b)
      return;

    CAPTURE(a);
    CAPTURE(b);
    CAPTURE(a - b);

    CAPTURE(diff_precision);
    if(diff_precision == 0)
      FAIL("diff_precision not initialized. "
           "Set diff_precision=-1 to enforce exact comparison, "
           "or to some positive number of bits.");

    if(diff_precision < 0)
      {
        REQUIRE(a == b);
        return;
      }

    REQUIRE(diff_precision > 0);

    auto eps = El::BigFloat(1) >>= diff_precision; // 2^{-precision}
    CAPTURE(eps);
    REQUIRE(Abs(a - b) < eps * (Abs(a) + Abs(b)));
  }
  template <class T1, class T2>
  void diff(const std::pair<T1, T2> &a, const std::pair<T1, T2> &b)
  {
    INFO("diff std::pair");
    DIFF(a.first, b.first);
    DIFF(a.second, b.second);
  }
  template <class T>
  void diff(const std::vector<T> &a, const std::vector<T> &b)
  {
    INFO("diff std::vector");
    REQUIRE(a.size() == b.size());
    for(size_t i = 0; i < a.size(); ++i)
      {
        CAPTURE(i);
        DIFF(a[i], b[i]);
      }
  }
  template <class T> void diff(const El::Matrix<T> &a, const El::Matrix<T> &b)
  {
    INFO("diff El::Matrix");
    REQUIRE(a.Height() == b.Height());
    REQUIRE(a.Width() == b.Width());
    for(int row = 0; row < a.Height(); ++row)
      for(int col = 0; col < a.Width(); ++col)
        {
          CAPTURE(row);
          CAPTURE(col);
          DIFF(a.Get(row, col), b.Get(row, col));
        }
  }
  template <class T>
  void diff(const El::DistMatrix<T> &a, const El::DistMatrix<T> &b,
            const std::optional<El::UpperOrLower> uplo = std::nullopt)
  {
    INFO("diff El::DistMatrix");
    CAPTURE(uplo);
    REQUIRE(a.Height() == b.Height());
    REQUIRE(a.Width() == b.Width());
    REQUIRE(a.Grid() == b.Grid());

    // only square matrices allowed for uplo
    if(uplo.has_value())
      REQUIRE(a.Height() == a.Width());

    for(int iLoc = 0; iLoc < a.LocalHeight(); ++iLoc)
      for(int jLoc = 0; jLoc < a.LocalWidth(); ++jLoc)
        {
          CAPTURE(iLoc);
          CAPTURE(jLoc);
          int i = a.GlobalRow(iLoc);
          int j = a.GlobalCol(jLoc);
          CAPTURE(i);
          CAPTURE(j);
          if(uplo.has_value())
            {
              if(uplo.value() == El::UPPER && i > j)
                continue;
              if(uplo.value() == El::LOWER && i < j)
                continue;
            }

          DIFF(a.GetLocal(iLoc, jLoc), b.GetLocal(iLoc, jLoc));
        }
  }
  template <>
  inline void
  diff(const std::filesystem::path &a, const std::filesystem::path &b)
  {
    CAPTURE(a);
    CAPTURE(b);
    DIFF(weakly_canonical(a).string(), weakly_canonical(b).string());
  }
  template <>
  inline void diff(const Damped_Rational &a, const Damped_Rational &b)
  {
    INFO("diff Damped_Rational");
    DIFF(a.base, b.base);
    DIFF(a.constant, b.constant);
    DIFF(a.poles, b.poles);
  }
  template <> inline void diff(const PVM_Info &a, const PVM_Info &b)
  {
    // Ignore path
    // DIFF(a.block_path,b.block_path);
    DIFF(a.prefactor, b.prefactor);
    DIFF(a.reduced_prefactor, b.reduced_prefactor);
    DIFF(a.sample_points, b.sample_points);
    DIFF(a.sample_scalings, b.sample_scalings);
    DIFF(a.reduced_sample_scalings, b.reduced_sample_scalings);
  }
}
