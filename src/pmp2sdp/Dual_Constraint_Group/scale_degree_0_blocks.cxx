// Two-set encoding (pmp2sdp --encodingPoints): the scale of the blocks of
// degree 0.
//
// A block of degree 0 (one sample point, e.g. an isolated operator) has
// nothing to re-encode: its rows are the block's values at its one point,
// at whatever scale the input gives them.  The re-encoded blocks' rows, on
// the other hand, are at the scale of their encoding scalings, which can be
// many orders of magnitude below the default encoding's (3d Ising with the
// stress tensor at nmax 22: ~50 orders).  A degree-0 block left at the
// input's scale then dominates SDPB's Schur complement Q at the initial
// point, and Cholesky(Q) fails in the first iteration at any precision the
// rest of the SDP needs.
//
// At the initial point X = Y = Omega I, a block's term in Q is B^T S^-1 B
// with S built from its bilinear bases.  The re-encoded blocks' bases are
// orthonormal (S of order 1), so their term is set by |B|.  A degree-0
// block's basis is the single number b = sqrt(sigma), sigma its sample
// scaling (A_p = b^2 E_rs, S = b^4), so its term is set by |B| / b^2.
//
// Multiplying a block's rows of B and c by a positive factor s, with its
// bilinear basis unchanged, gives an equivalent constraint:
// Tr(A_p Y) + s (B y)_p = s c_p is Tr(A_p Y / s) + (B y)_p = c_p.  The
// feasible set in y, the dual and primal objectives and the optimum are
// unchanged; only that block's X and x are divided, and its Y multiplied, by
// s.  (Multiplying sigma by s instead would scale A_p with B and c and leave
// the block's term in Q as it was.)  So pmp2sdp puts every degree-0 block at
// the scale of the re-encoded blocks: s = 2^e, with e the integer nearest to
// log2(M / m), where m is the block's largest |entry| of B and c divided by
// b^2, and M the median over the re-encoded blocks (all ranks) of their
// largest |entry| of B and c.  A power of 2 changes no digit of the rows.
// The block's sample scalings (written to pmp_info.json) are multiplied by s
// too, so that its rows remain "scaling times polynomial" and x times the
// scaling is unchanged, which is what spectrum reads back.
//
// Why the median: SDPB starts every block at the same X and Y (multiples of
// the identity), so a block whose rows are far from all the others' -- above
// or below -- costs iterations, and the re-encoded blocks' own scales can
// span tens of orders of magnitude.  The median is the scale least far from
// all of them and is not moved by a few extreme blocks.  (Measured on a
// scalar nmax 14 gap problem with one constant block: two-set needs the same
// number of iterations as without the block when the block is anywhere inside
// the range of the re-encoded blocks' scales, and more when it is above or
// below it; at the input's scale Cholesky(Q) fails.)

#include "../Dual_Constraint_Group.hxx"
#include "sdpb_util/assert.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

namespace
{
  // log2 |x| (x != 0)
  double log2_abs(const El::BigFloat &x)
  {
    long exponent;
    const double mantissa = mpf_get_d_2exp(&exponent, x.gmp_float.get_mpf_t());
    return static_cast<double>(exponent) + std::log2(std::abs(mantissa));
  }

  // log2 of the largest |entry| of B and c; none if all are zero
  std::optional<double> log2_max_abs(const Dual_Constraint_Group &group)
  {
    std::optional<double> result;
    auto update = [&result](const El::BigFloat &x) {
      if(mpf_sgn(x.gmp_float.get_mpf_t()) == 0)
        return;
      const double value = log2_abs(x);
      if(!result.has_value() || value > result.value())
        result = value;
    };
    for(const auto &c : group.constraint_constants)
      update(c);
    const auto &B = group.constraint_matrix;
    for(El::Int j = 0; j < B.Width(); ++j)
      for(El::Int i = 0; i < B.Height(); ++i)
        update(B.CRef(i, j));
    return result;
  }

  void multiply_by_power_of_2(El::BigFloat &x, const int64_t exponent)
  {
    auto *const f = x.gmp_float.get_mpf_t();
    if(exponent > 0)
      mpf_mul_2exp(f, f, static_cast<mp_bitcnt_t>(exponent));
    else if(exponent < 0)
      mpf_div_2exp(f, f, static_cast<mp_bitcnt_t>(-exponent));
  }

  // All ranks' values, in rank order
  std::vector<double>
  all_gather(const std::vector<double> &local, const El::mpi::Comm &comm)
  {
    const int num_procs = El::mpi::Size(comm);
    const int local_count = static_cast<int>(local.size());
    std::vector<int> counts(num_procs), displacements(num_procs, 0);
    MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
                  comm.comm);
    for(int rank = 1; rank < num_procs; ++rank)
      displacements[rank] = displacements[rank - 1] + counts[rank - 1];
    std::vector<double> result(displacements.back() + counts.back());
    MPI_Allgatherv(local.data(), local_count, MPI_DOUBLE, result.data(),
                   counts.data(), displacements.data(), MPI_DOUBLE,
                   comm.comm);
    return result;
  }
}

double median_log2_max_abs(std::vector<double> values)
{
  ASSERT(!values.empty());
  std::sort(values.begin(), values.end());
  const size_t n = values.size();
  return n % 2 == 1 ? values[n / 2] : (values[n / 2 - 1] + values[n / 2]) / 2;
}

Degree_0_Scaling
scale_degree_0_blocks(std::vector<Dual_Constraint_Group> &groups,
                      std::vector<Polynomial_Vector_Matrix> &matrices,
                      const El::mpi::Comm &comm)
{
  ASSERT_EQUAL(groups.size(), matrices.size());
  std::vector<double> reencoded;
  std::vector<std::optional<double>> scales(groups.size());
  El::Int num_degree_0 = 0;
  for(size_t i = 0; i < groups.size(); ++i)
    {
      ASSERT_EQUAL(groups[i].num_points, matrices[i].sample_points.size());
      scales[i] = log2_max_abs(groups[i]);
      if(groups[i].num_points < 2)
        {
          ++num_degree_0;
          // relative to the block's A = b^2 (see above)
          const auto &basis = groups[i].bilinear_bases[0];
          ASSERT(basis.Height() == 1 && basis.Width() == 1,
                 "a degree-0 block has a 1x1 even bilinear basis",
                 DEBUG_STRING(basis.Height()), DEBUG_STRING(basis.Width()));
          ASSERT(mpf_sgn(basis.CRef(0, 0).gmp_float.get_mpf_t()) != 0);
          if(scales[i].has_value())
            scales[i] = scales[i].value() - 2 * log2_abs(basis.CRef(0, 0));
        }
      else if(scales[i].has_value())
        reencoded.push_back(scales[i].value());
    }

  Degree_0_Scaling result;
  result.num_degree_0 = El::mpi::AllReduce(num_degree_0, El::mpi::SUM, comm);
  const auto all_reencoded = all_gather(reencoded, comm);
  result.num_reencoded = static_cast<El::Int>(all_reencoded.size());
  if(result.num_degree_0 == 0 || all_reencoded.empty())
    return result;
  result.median_log2 = median_log2_max_abs(all_reencoded);

  El::Int min_exponent = std::numeric_limits<El::Int>::max(),
          max_exponent = std::numeric_limits<El::Int>::min(), num_scaled = 0;
  for(size_t i = 0; i < groups.size(); ++i)
    {
      if(groups[i].num_points >= 2 || !scales[i].has_value())
        continue;
      const El::Int exponent
        = std::llround(result.median_log2 - scales[i].value());
      min_exponent = std::min(min_exponent, exponent);
      max_exponent = std::max(max_exponent, exponent);
      if(exponent == 0)
        continue;
      ++num_scaled;
      auto &group = groups[i];
      for(auto &c : group.constraint_constants)
        multiply_by_power_of_2(c, exponent);
      for(El::Int j = 0; j < group.constraint_matrix.Width(); ++j)
        for(El::Int k = 0; k < group.constraint_matrix.Height(); ++k)
          multiply_by_power_of_2(group.constraint_matrix.Ref(k, j), exponent);
      auto &m = matrices[i];
      for(auto &s : m.sample_scalings)
        multiply_by_power_of_2(s, exponent);
      for(auto &s : m.reduced_sample_scalings)
        multiply_by_power_of_2(s, exponent);
    }
  result.num_scaled = El::mpi::AllReduce(num_scaled, El::mpi::SUM, comm);
  result.min_exponent = El::mpi::AllReduce(min_exponent, El::mpi::MIN, comm);
  result.max_exponent = El::mpi::AllReduce(max_exponent, El::mpi::MAX, comm);
  if(result.min_exponent > result.max_exponent)
    result.min_exponent = result.max_exponent = 0; // only zero blocks
  return result;
}
