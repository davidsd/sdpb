#include <catch2/catch_amalgamated.hpp>
#include <boost_serialization.hxx>
#include <El.hpp>
#include <sstream>
#include "unit_tests/util/util.hxx"

using Test_Util::random_bigfloat;
using Test_Util::random_matrix;
using Test_Util::REQUIRE_Equal::diff;

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
    DIFF(matrix, other);
  }
}
