#pragma once

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
    inline explicit Diff_Precision(int precision)
    {
      old_precision = diff_precision;
      diff_precision = precision;
    }
    inline virtual ~Diff_Precision() { diff_precision = old_precision; }

  private:
    int old_precision;
  };

  inline void diff(int a, int b)
  {
    REQUIRE(a == b);
  }
  inline void diff(const El::BigFloat &a, const El::BigFloat &b)
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
  inline void diff(const std::string &a, const std::string &b)
  {
    REQUIRE(a == b);
  }
  template <class T1, class T2>
  inline void diff(const std::pair<T1, T2> &a, const std::pair<T1, T2> &b)
  {
    INFO("diff std::pair");
    DIFF(a.first, b.first);
    DIFF(a.second, b.second);
  }
  template <class T>
  inline void diff(const std::vector<T> &a, const std::vector<T> &b)
  {
    INFO("diff std::vector");
    REQUIRE(a.size() == b.size());
    for(size_t i = 0; i < a.size(); ++i)
      {
        CAPTURE(i);
        DIFF(a[i], b[i]);
      }
  }
  template <class T>
  inline void diff(const El::Matrix<T> &a, const El::Matrix<T> &b)
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
  inline void diff(const El::DistMatrix<T> &a, const El::DistMatrix<T> &b)
  {
    INFO("diff El::DistMatrix");
    REQUIRE(a.Height() == b.Height());
    REQUIRE(a.Width() == b.Width());
    REQUIRE(a.Grid() == b.Grid());

    for(int row = 0; row < a.LocalHeight(); ++row)
      for(int col = 0; col < a.LocalWidth(); ++col)
        {
          CAPTURE(row);
          CAPTURE(col);
          DIFF(a.GetLocal(row, col), b.GetLocal(row, col));
        }
  }
}
