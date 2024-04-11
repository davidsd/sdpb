#pragma once

#include "test_util/test_util.hxx"

#include <vector>
#include <El.hpp>

namespace Test_Util
{
  inline El::BigFloat random_bigfloat()
  {
    return El::SampleUniform<El::BigFloat>(-1.0, 1.0);
  }

  inline std::vector<El::BigFloat>
  random_vector(size_t size, const std::function<El::BigFloat()> &make_float
                             = random_bigfloat)
  {
    std::vector<El::BigFloat> vec(size);
    for(auto &item : vec)
      {
        item = make_float();
      }
    return vec;
  }

  inline El::Matrix<El::BigFloat>
  random_matrix(int height, int width,
                const std::function<El::BigFloat()> &make_float
                = random_bigfloat)
  {
    El::Matrix<El::BigFloat> matrix(height, width);
    for(int i = 0; i < height; ++i)
      for(int k = 0; k < width; ++k)
        {
          matrix.Set(i, k, make_float());
        }

    return matrix;
  }

  inline El::Matrix<El::BigFloat> zero_matrix(int height, int width)
  {
    El::Matrix<El::BigFloat> zeros(height, width);
    El::Zero(zeros);
    return zeros;
  }

  inline El::DistMatrix<El::BigFloat>
  random_distmatrix(int height, int width,
                    const std::function<El::BigFloat()> &make_float
                    = random_bigfloat)
  {
    El::DistMatrix<El::BigFloat> matrix(height, width);
    for(int i = 0; i < matrix.LocalHeight(); ++i)
      for(int k = 0; k < matrix.LocalWidth(); ++k)
        {
          matrix.SetLocal(i, k, make_float());
        }

    return matrix;
  }

  // Split n>0 into sum of random positive integers n1+n2+...
  template <class TInt> inline std::vector<TInt> random_split(TInt n)
  {
    std::vector<TInt> result;

    std::default_random_engine rand_engine;

    while(n > 0)
      {
        std::uniform_int_distribution<TInt> dist(1, n);
        int curr = dist(rand_engine);
        result.push_back(curr);
        n -= curr;
      }
    return result;
  }
}
