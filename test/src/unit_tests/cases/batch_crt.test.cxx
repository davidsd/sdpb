#include "catch2/catch_amalgamated.hpp"

#include "sdpb_util/bigint_shared_memory/fmpz/Batch_Crt.hxx"
#include "sdpb_util/bigint_shared_memory/fmpz/Fmpz_BigInt.hxx"
#include "sdpb_util/bigint_shared_memory/fmpz/residues.hxx"
#include "unit_tests/util/util.hxx"

#include <El.hpp>

#include <cmath>
#include <cstring>
#include <random>
#include <vector>

// Batch_Crt (sdpb_util/bigint_shared_memory/fmpz/Batch_Crt.hxx) against the
// per-entry path of restore_and_reduce (double_to_uint32_t_residue,
// fmpz_multi_CRT_ui, fmpz_get_mpf, El::BigFloat::Serialize), byte for byte in
// the serialized size, exponent and limbs.

namespace
{
  void random_fmpz(fmpz_t x, std::mt19937_64 &rng, const long bits)
  {
    fmpz_zero(x);
    for(long b = 0; b < bits; b += 64)
      {
        fmpz_mul_2exp(x, x, 64);
        fmpz_add_ui(x, x, rng());
      }
    fmpz_fdiv_r_2exp(x, x, bits);
    if(rng() & 1)
      fmpz_neg(x, x);
  }

  // Meaningful bytes of a serialized BigFloat: size, exponent, |size| limbs
  std::vector<El::byte> meaningful(const El::byte *slot)
  {
    int size;
    std::memcpy(&size, slot, sizeof(int));
    const size_t n
      = sizeof(int) + sizeof(mp_exp_t) + std::abs(size) * sizeof(mp_limb_t);
    return std::vector<El::byte>(slot, slot + n);
  }
}

TEST_CASE("batch_crt")
{
  if(El::mpi::Rank() != 0)
    return;
  std::mt19937_64 rng(269);
  // bits of the comb's factors (the syrk uses El::gmp::Precision()); the
  // BigFloats are at the working precision, so integers longer than it are
  // truncated as mpf_set_z does
  const long working = El::gmp::Precision();
  for(const long factor_bits : {64L, 256L, working, 1024L})
    for(const slong k : {slong(1), slong(13), slong(5490)})
      for(const double band : {1e-9, 1.0})
        DYNAMIC_SECTION("factor_bits=" << factor_bits << " k=" << k
                                       << " fallback_band=" << band)
        {
          // as BigInt_Shared_Memory_Syrk_Context: comb(prec, prec, 1, k)
          Fmpz_Comb comb(factor_bits, factor_bits, 1, k);
          // products below 200 bits with small k: FLINT's primes are too
          // large for the exact sums; restore_and_reduce then keeps the
          // per-entry path
          if(factor_bits == 64 && k == 1)
            REQUIRE(!Batch_Crt::supported(comb));
          if(factor_bits >= 256)
            REQUIRE(Batch_Crt::supported(comb));
          if(!Batch_Crt::supported(comb))
            REQUIRE_THROWS(Batch_Crt(comb));
          else
          {
          Batch_Crt crt(comb);
          crt.fallback_band = band;
          const size_t K = comb.num_primes;
          // |c| < k 2^(2 factor_bits) < M/2
          const long max_bits = 2 * factor_bits + FLINT_BIT_COUNT(k) - 1;

          // 0; +-2^b and +-(2^b - 1) for b = 0, 37, 74, ...; random sizes
          std::vector<fmpz> c;
          c.push_back(0);
          for(long b = 0; b < max_bits; b += 37)
            for(int s = 0; s < 4; ++s)
              {
                fmpz_t v;
                fmpz_init(v);
                fmpz_one(v);
                fmpz_mul_2exp(v, v, b);
                if(s % 2 == 1)
                  fmpz_sub_ui(v, v, 1);
                if(s >= 2)
                  fmpz_neg(v, v);
                c.push_back(*v);
              }
          for(int r = 0; r < 2000; ++r)
            {
              fmpz_t v;
              fmpz_init(v);
              random_fmpz(
                v, rng, std::uniform_int_distribution<long>(1, max_bits)(rng));
              c.push_back(*v);
            }
          const size_t n = c.size();
          // residues as the BLAS jobs leave them: any representative
          // r = (c mod p) + m p, |r| <= 2^53
          std::vector<std::vector<double>> residues(K, std::vector<double>(n));
          for(size_t i = 0; i < K; ++i)
            {
              const mp_limb_t p = comb.primes[i];
              const double m_max = std::floor(std::ldexp(1.0, 53) / p) - 1;
              std::uniform_real_distribution<double> m(-m_max, m_max);
              for(size_t e = 0; e < n; ++e)
                residues[i][e] = double(fmpz_fdiv_ui(&c[e], p))
                                 + std::floor(m(rng)) * double(p);
            }
          std::vector<const double *> ptr(K);
          for(size_t i = 0; i < K; ++i)
            ptr[i] = residues[i].data();
          const El::BigFloat zero;
          const size_t slot = zero.SerializedSize();
          const int mpf_prec = zero.gmp_float.get_mpf_t()[0]._mp_prec;
          std::vector<El::byte> out(n * slot);
          const size_t fallbacks = crt.restore_serialized(
            ptr.data(), n, comb, mpf_prec, out.data(), slot);
          if(band >= 1.0)
            REQUIRE(fallbacks == n);
          else
            CHECK(fallbacks == 0);

          // reference: residues_to_bigint's CRT + to_BigFloat + Serialize
          Fmpz_BigInt reference;
          El::BigFloat value;
          std::vector<El::byte> buf(slot);
          std::vector<mp_limb_t> res(K);
          for(size_t e = 0; e < n; ++e)
            {
              for(size_t i = 0; i < K; ++i)
                res[i]
                  = double_to_uint32_t_residue(residues[i][e], comb.mods[i]);
              fmpz_multi_CRT_ui(reference.value, res.data(), comb.comb,
                                comb.comb_temp, 1);
              CAPTURE(e);
              REQUIRE(fmpz_equal(reference.value, &c[e]));
              reference.to_BigFloat(value);
              value.Serialize(buf.data());
              REQUIRE(meaningful(buf.data())
                      == meaningful(out.data() + e * slot));
            }

          // negative control: one residue changed by 1 changes the result
          const size_t e = n - 1;
          std::vector<double> changed(K);
          std::vector<const double *> ptr1(K);
          for(size_t i = 0; i < K; ++i)
            {
              changed[i] = residues[i][e] + (i == 0 ? 1 : 0);
              ptr1[i] = &changed[i];
            }
          std::vector<El::byte> out1(slot);
          crt.restore_serialized(ptr1.data(), 1, comb, mpf_prec, out1.data(),
                                 slot);
          REQUIRE(meaningful(out1.data()) != meaningful(out.data() + e * slot));

          for(auto &v : c)
            fmpz_clear(&v);
          }
        }
}

// Batch_Crt::residues against Fmpz_BigInt::from_BigFloat +
// bigint_to_residues (compute_block_residues' per-entry path):
// the same doubles, for normalized-and-shifted entries |x| <= 2^P with
// fractional parts, of both signs, zero and exactly +-2^P.
TEST_CASE("batch_residues")
{
  if(El::mpi::Rank() != 0)
    return;
  std::mt19937_64 rng(2690);
  const long P = El::gmp::Precision();
  for(const slong k : {slong(13), slong(5490), slong(20000)})
    DYNAMIC_SECTION("k=" << k)
    {
      Fmpz_Comb comb(P, P, 1, k);
      REQUIRE(Batch_Crt::supported(comb));
      Batch_Crt crt(comb);
      const size_t K = comb.num_primes;
      std::vector<El::BigFloat> v;
      v.emplace_back(0);
      El::BigFloat two_P(1);
      two_P <<= P;
      v.push_back(two_P);
      v.push_back(-two_P);
      v.emplace_back(1);
      v.emplace_back(-1);
      v.emplace_back(0.5);
      v.emplace_back(-0.75);
      for(int r = 0; r < 3000; ++r)
        {
          // random in (-1, 1), times 2^(P - s), s uniform in [0, P + 40]
          El::BigFloat x = Test_Util::random_bigfloat();
          const long sh = P - std::uniform_int_distribution<long>(0, P + 40)(rng);
          if(sh >= 0)
            x <<= sh;
          else
            x >>= -sh;
          v.push_back(x);
        }
      const size_t n = v.size();
      std::vector<const __mpf_struct *> x(n);
      for(size_t e = 0; e < n; ++e)
        x[e] = v[e].gmp_float.get_mpf_t();
      // with a row stride, as the elementwise path's scatter would see it
      for(const size_t row_stride : {size_t(1), size_t(3)})
        {
          CAPTURE(row_stride);
          const size_t prime_stride = n * row_stride + 1;
          std::vector<double> out(prime_stride * K, -7.0);
          crt.residues(x.data(), n, P + 2, out.data(), row_stride,
                       prime_stride);
          Fmpz_BigInt bigint;
          std::vector<double> ref(K);
          for(size_t e = 0; e < n; ++e)
            {
              CAPTURE(e, v[e]);
              bigint.from_BigFloat(v[e]);
              bigint_to_residues(bigint, comb, ref.data(), 1);
              for(size_t i = 0; i < K; ++i)
                {
                  CAPTURE(i);
                  REQUIRE(out[e * row_stride + i * prime_stride] == ref[i]);
                }
            }
        }
      // an input above 2^bits is refused
      El::BigFloat big(1);
      big <<= P + 3;
      const __mpf_struct *xb = big.gmp_float.get_mpf_t();
      std::vector<double> o(K);
      REQUIRE_THROWS(crt.residues(&xb, 1, P + 2, o.data(), 1, 1));
    }
}
