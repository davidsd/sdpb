#include <catch2/catch_amalgamated.hpp>
#include <boost_serialization.hxx>
#include <El.hpp>
#include <sstream>

namespace
{
  template <class T> T serialize_deserialize(const T &value)
  {
    std::stringstream ss;
    boost::archive::binary_oarchive out(ss);
    out << value;

    T deserialized_value;
    boost::archive::binary_iarchive in(ss);
    in >> deserialized_value;

    return deserialized_value;
  }

  El::BigFloat random_bigfloat()
  {
    return El::SampleUniform<El::BigFloat>(-3.14, 3.14);
  }
  El::Matrix<El::BigFloat> random_matrix(int height, int width)
  {
    El::Matrix<El::BigFloat> matrix(height, width);
    for(int i = 0; i < height; ++i)
      for(int k = 0; k < width; ++k)
        {
          matrix.Set(i, k, random_bigfloat());
        }

    return matrix;
  }

  void diff(El::Matrix<El::BigFloat> a, El::Matrix<El::BigFloat> b)
  {
    // Compare dimensions
    REQUIRE(b.Height() == a.Height());
    REQUIRE(a.Width() == b.Width());
    REQUIRE(a.LDim() == b.LDim());
    // Elementwise comparison
    for(int i = 0; i < a.Height(); ++i)
      for(int k = 0; k < a.Width(); ++k)
        REQUIRE(a.Get(i, k) == b.Get(i, k));
  }
}

TEST_CASE("Boost serialization")
{
  // this test is purely single-process, thus testing at rank=0 is sufficient
  if(El::mpi::Rank() != 0)
    return;

  El::InitializeRandom(true);

  SECTION("El::BigFloat")
  {
    El::BigFloat value = random_bigfloat();
    El::BigFloat other = serialize_deserialize(value);
    REQUIRE(value == other);
  }

  SECTION("El::Matrix<BigFloat>")
  {
    int height = 2;
    int width = 3;
    auto matrix = random_matrix(height, width);

    El::Matrix<El::BigFloat> other = serialize_deserialize(matrix);
    // Sanity check: deserialized_matrix is not the same as matrix
    REQUIRE(matrix.LockedBuffer() != other.LockedBuffer());
    diff(matrix, other);
  }
}
