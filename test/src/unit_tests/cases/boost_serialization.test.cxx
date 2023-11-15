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
    auto zero = El::BigFloat(0);
    auto nonzero = random_bigfloat();

    for(auto &value : {zero, nonzero})
      {
        CAPTURE(value);
        El::BigFloat other = serialize_deserialize(value);
        REQUIRE(value == other);
      }

    {
      INFO("Check that zero serialization is more compact");
      std::stringstream ss_zero, ss_nonzero;
      boost::archive::binary_oarchive ar_zero(ss_zero), ar_nonzero(ss_nonzero);
      ar_zero << zero;
      ar_nonzero << nonzero;
      REQUIRE(ss_zero.str().size() < ss_nonzero.str().size());
    }
  }

  SECTION("El::Matrix<BigFloat>")
  {
    int height = 100;
    int width = 10;
    auto rand_matrix = random_matrix(height, width);
    El::Matrix<El::BigFloat> zeros(height, width);
    El::Zero(zeros);

    for(auto &matrix : {zeros, rand_matrix})
      {
        El::Matrix<El::BigFloat> other = serialize_deserialize(matrix);
        // Sanity check: deserialized_matrix is not the same as matrix
        REQUIRE(matrix.LockedBuffer() != other.LockedBuffer());
        DIFF(matrix, other);
      }
  }
}
