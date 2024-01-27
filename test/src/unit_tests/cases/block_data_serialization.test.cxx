#include <catch2/catch_amalgamated.hpp>
#include <El.hpp>
#include <sstream>

#include "unit_tests/util/util.hxx"
#include "pmp2sdp/Dual_Constraint_Group.hxx"
#include "pmp2sdp/Block_File_Format.hxx"

using Test_Util::random_bigfloat;
using Test_Util::random_matrix;
using Test_Util::random_vector;
using Test_Util::zero_matrix;
using Test_Util::REQUIRE_Equal::diff;

void write_block_data(std::ostream &os, const Dual_Constraint_Group &group,
                      Block_File_Format format);
void parse_block_data(std::istream &block_stream, Block_File_Format format,
                      El::Matrix<El::BigFloat> &constraint_matrix,
                      std::vector<El::BigFloat> &constraint_constants,
                      El::Matrix<El::BigFloat> &bilinear_bases_even,
                      El::Matrix<El::BigFloat> &bilinear_bases_odd);

namespace
{
  Dual_Constraint_Group
  serialize_deserialize(const Dual_Constraint_Group &group,
                        Block_File_Format format)
  {
    std::stringstream ss;
    write_block_data(ss, group, format);

    Dual_Constraint_Group result;
    parse_block_data(ss, format, result.constraint_matrix,
                     result.constraint_constants, result.bilinear_bases[0],
                     result.bilinear_bases[1]);
    result.dim = group.dim;
    result.num_points = group.num_points;
    return result;
  }

  using Test_Util::REQUIRE_Equal::diff;
  void diff(const Dual_Constraint_Group &a, const Dual_Constraint_Group &b)
  {
    DIFF(a.dim, b.dim);
    DIFF(a.num_points, b.num_points);
    DIFF(a.constraint_matrix, b.constraint_matrix);
    DIFF(a.constraint_constants, b.constraint_constants);
    DIFF(a.bilinear_bases[0], b.bilinear_bases[0]);
    DIFF(a.bilinear_bases[1], b.bilinear_bases[1]);
  }

  // Dual_Constraint_Group with the same sizes as block_0 in integration test
  // end-to-end_tests/SingletScalar_cT_test_nmax6/primal_dual_optimal
  Dual_Constraint_Group random_group_from_singlet_scalar_block_0()
  {
    Dual_Constraint_Group group;
    group.dim = 1;
    group.num_points = 24;
    int P = 24;
    int N = 20;
    group.constraint_matrix = random_matrix(P, N);
    group.constraint_constants = random_vector(P);
    group.bilinear_bases[0] = random_matrix(12, 24);
    group.bilinear_bases[1] = random_matrix(12, 24);
    return group;
  }

  // Dual_Constraint_Group with the same sizes as block_0 in integration test
  // end-to-end_tests/SingletScalar_cT_test_nmax6/primal_dual_optimal
  Dual_Constraint_Group zero_group_from_singlet_scalar_block_0()
  {
    Dual_Constraint_Group group;
    group.dim = 1;
    group.num_points = 24;
    int P = 24;
    int N = 20;
    group.constraint_matrix = zero_matrix(P, N);
    group.constraint_constants = std::vector<El::BigFloat>(P, 0);
    group.bilinear_bases[0] = zero_matrix(12, 24);
    group.bilinear_bases[1] = zero_matrix(12, 24);
    return group;
  }
}

TEST_CASE("block_data serialization")
{
  // this test is purely single-process, thus testing at rank=0 is sufficient
  if(El::mpi::Rank() != 0)
    return;

  El::InitializeRandom(true);
  Dual_Constraint_Group group = random_group_from_singlet_scalar_block_0();
  Dual_Constraint_Group zero_group = zero_group_from_singlet_scalar_block_0();

  Block_File_Format format = GENERATE(bin, json);
  DYNAMIC_SECTION((format == bin ? ".bin" : ".json"))
  {
    auto other = serialize_deserialize(group, format);
    DIFF(group, other);
    other = serialize_deserialize(zero_group, format);
    DIFF(zero_group, other);
  }
}

TEST_CASE("benchmark block_data write+parse", "[!benchmark]")
{
  // this test is purely single-process, thus testing at rank=0 is sufficient
  if(El::mpi::Rank() != 0)
    return;

  El::InitializeRandom(true);
  Dual_Constraint_Group group = random_group_from_singlet_scalar_block_0();
  Dual_Constraint_Group zero_group = zero_group_from_singlet_scalar_block_0();

  // Change constraint_matrix size to see how bin/json scales
  int B_width = GENERATE(20, 100, 1000, 10000);
  group.constraint_matrix
    = random_matrix(group.constraint_matrix.Height(), B_width);
  zero_group.constraint_matrix
    = zero_matrix(zero_group.constraint_matrix.Height(), B_width);

  int total_count = 0;
  total_count += group.constraint_constants.size();
  for(const auto &matrix : {group.constraint_matrix, group.bilinear_bases[0],
                            group.bilinear_bases[1]})
    {
      total_count += matrix.Height() * matrix.Width();
    }

  DYNAMIC_SECTION(total_count << " BigFloats")
  {
    // We could put this benchmarks into different DYNAMIC_SECTION's,
    // using Block_File_Format format = GENERATE(bin, json);
    // But it would make output less concise
    BENCHMARK("write+parse bin zero")
    {
      return serialize_deserialize(zero_group, bin);
    };
    BENCHMARK("write+parse bin nonzero")
    {
      return serialize_deserialize(group, bin);
    };
    BENCHMARK("write+parse JSON zero")
    {
      return serialize_deserialize(zero_group, json);
    };
    BENCHMARK("write+parse JSON nonzero")
    {
      return serialize_deserialize(group, json);
    };
  }
}
