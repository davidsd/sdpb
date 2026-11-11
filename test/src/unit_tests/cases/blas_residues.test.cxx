#include "catch2/catch_amalgamated.hpp"

#include "sdpb_util/bigint_shared_memory/fmpz/blas_residues.hxx"
#include "sdpb_util/bigint_shared_memory/fmpz/residues.hxx"
#include "unit_tests/util/util.hxx"

#include <El.hpp>

#include <climits>
#include <random>
#include <vector>

// BLAS residues and BLAS CRT (sdpb_util/bigint_shared_memory/fmpz/blas_residues.hxx) against
// FLINT's per-element routines.

namespace
{
  // Random BigFloat in (-1, 1) times 2^e, e uniform in [-max_shift, 0]
  // (entries of very different sizes in one matrix)
  El::BigFloat random_scaled(std::mt19937_64 &rng, const int max_shift)
  {
    El::BigFloat x = Test_Util::random_bigfloat();
    const int shift = std::uniform_int_distribution<int>(0, max_shift)(rng);
    mpf_div_2exp(x.gmp_float.get_mpf_t(), x.gmp_float.get_mpf_t(), shift);
    return x;
  }

  // Random integer of at most `bits` bits, random sign
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

  // Integer trunc(x 2^shift) as an fmpz
  void fixed_point_fmpz(const El::BigFloat &x, const long shift, fmpz_t out)
  {
    mpf_t scaled;
    mpf_init2(scaled, mpf_get_prec(x.gmp_float.get_mpf_t()) + 64);
    if(shift >= 0)
      mpf_mul_2exp(scaled, x.gmp_float.get_mpf_t(), shift);
    else
      mpf_div_2exp(scaled, x.gmp_float.get_mpf_t(), -shift);
    mpz_t z;
    mpz_init(z);
    mpz_set_f(z, scaled); // truncates toward zero
    fmpz_set_mpz(out, z);
    mpz_clear(z);
    mpf_clear(scaled);
  }
}

TEST_CASE("blas_residues")
{
  if(El::mpi::Rank() != 0)
    return;
  std::mt19937_64 rng(1);
  // bits < 141: exact BLAS CRT (t0 = 0); bits = precision + 64 guard bits:
  // the setting of the multimodular Cholesky/Trsm
  const long bits
    = GENERATE(long(110), long(140), long(El::gmp::Precision() + 64));
  const slong k = GENERATE(slong(1), slong(37), slong(1000));
  DYNAMIC_SECTION("bits=" << bits << " k=" << k)
  {
    const auto &context = Blas_Residues::get_context(bits, bits, k);
    const size_t num_primes = context.num_primes;
    CAPTURE(num_primes);
    CAPTURE(context.t0);
    if(bits < 141)
      REQUIRE(context.t0 == 0);

    SECTION("residues equal FLINT's")
    {
      // Fixed-point integers trunc(x 2^(bits - e_max)) of random x of
      // different sizes, so that |x 2^shift| < 2^bits
      const size_t n = 3000; // more than one batch of 2048
      std::vector<El::BigFloat> values(n);
      for(auto &value : values)
        value = random_scaled(rng, 200);
      values[5] = 0;
      long e_max = LONG_MIN;
      for(const auto &value : values)
        if(value != El::BigFloat(0))
          {
            long e;
            mpf_get_d_2exp(&e, value.gmp_float.get_mpf_t());
            e_max = std::max(e_max, e);
          }
      const long shift = bits - e_max;

      std::vector<int32_t> residues(n * num_primes);
      for(const bool multiply_by_u : {false, true})
        {
          CAPTURE(multiply_by_u);
          Blas_Residues::compute_blas_residues(
            context, n,
            [&](size_t e, long &s) -> const __mpf_struct * {
              s = shift;
              return values[e].gmp_float.get_mpf_t();
            },
            multiply_by_u, residues.data());

          Fmpz_BigInt x;
          std::vector<double> flint_residues(num_primes);
          for(size_t e = 0; e < n; ++e)
            {
              fixed_point_fmpz(values[e], shift, x.value);
              REQUIRE(fmpz_bits(x.value) <= flint_bitcnt_t(bits));
              // FLINT >= 3.7 takes a non-const comb; it writes only to its scratch space
              bigint_to_residues(x, const_cast<Fmpz_Comb &>(context.comb),
                                 flint_residues.data(), 1);
              for(size_t p = 0; p < num_primes; ++p)
                {
                  const auto prime = context.comb.primes[p];
                  // FLINT's residue in (-p/2, p/2], times u_p if required
                  int64_t expected = int64_t(flint_residues[p]);
                  if(multiply_by_u)
                    expected
                      = (expected * int64_t(context.u[p])) % int64_t(prime);
                  const int64_t r = residues[e + n * p];
                  CAPTURE(e);
                  CAPTURE(p);
                  REQUIRE(std::abs(r) <= int64_t(prime + 1) / 2);
                  REQUIRE(((r - expected) % int64_t(prime)) == 0);
                }
            }
        }
    }

    SECTION("BLAS CRT vs FLINT CRT")
    {
      // Random products c, |c| < k 2^(2 bits), as from multimodular GEMM
      const size_t n = 1500; // more than one batch
      const long c_bits = 2 * bits + FLINT_BIT_COUNT(k);
      std::vector<fmpz> c(n);
      for(size_t e = 0; e < n; ++e)
        {
          fmpz_init(&c[e]);
          // a few small values and zero, the rest of up to full size
          const long b = e % 7 == 0 ? (long)(e % 50) : c_bits - 1;
          random_fmpz(&c[e], rng, b);
        }
      fmpz_zero(&c[3]);

      // y_p = c u_p mod p in [0, p) for the BLAS CRT, c mod p for FLINT,
      // and quotient = sum_p y_p / p as accumulated in multimodular GEMM
      std::vector<int32_t> y(n * num_primes), y_flint(n * num_primes);
      std::vector<double> quotients(n, 0.0);
      for(size_t p = 0; p < num_primes; ++p)
        {
          const auto prime = context.comb.primes[p];
          for(size_t e = 0; e < n; ++e)
            {
              const ulong r = fmpz_fdiv_ui(&c[e], prime);
              y_flint[e + n * p] = r;
              const double r_u = double((r * ulong(context.u[p])) % prime);
              y[e + n * p] = int32_t(r_u);
              quotients[e] += r_u * context.inverse_primes[p];
            }
        }

      Blas_Residues::Crt_Workspace work(context);
      Blas_Residues::blas_crt_high_chunks(context, y.data(), n, n, work);
      mpf_t blas_value, flint_value, diff, bound_low, bound_high;
      const mp_bitcnt_t exact_precision = c_bits + 128;
      for(auto *v :
          {&blas_value, &flint_value, &diff, &bound_low, &bound_high})
        mpf_init2(*v, exact_precision);
      // c - c' must lie in (-2^(24 t0), num_primes 2^21 2^(24 t0))
      mpf_set_ui(bound_low, 1);
      mpf_mul_2exp(bound_low, bound_low,
                   Blas_Residues::chunk_bits * context.t0);
      mpf_neg(bound_low, bound_low);
      mpf_set_ui(bound_high, num_primes);
      mpf_mul_2exp(bound_high, bound_high,
                   21 + Blas_Residues::chunk_bits * context.t0);
      for(size_t e = 0; e < n; ++e)
        {
          CAPTURE(e);
          mpf_set_ui(blas_value, 0);
          mpf_set_ui(flint_value, 0);
          Blas_Residues::blas_crt_to_mpf(context, work.S.data() + e, n,
                                         quotients[e], 0, work, blas_value);
          Blas_Residues::flint_crt_to_mpf(context, y_flint.data() + e, n, 0,
                                          work, flint_value);
          // FLINT's CRT is exact
          mpz_t z;
          mpz_init(z);
          fmpz_get_mpz(z, &c[e]);
          mpf_set_z(diff, z);
          mpz_clear(z);
          REQUIRE(mpf_cmp(flint_value, diff) == 0);

          mpf_sub(diff, flint_value, blas_value);
          if(context.t0 == 0)
            REQUIRE(mpf_sgn(diff) == 0);
          else
            {
              REQUIRE(mpf_cmp(diff, bound_low) > 0);
              REQUIRE(mpf_cmp(diff, bound_high) < 0);
            }
        }
      for(auto *v :
          {&blas_value, &flint_value, &diff, &bound_low, &bound_high})
        mpf_clear(*v);
      for(auto &x : c)
        fmpz_clear(&x);
    }
  }
}

// Exact dyadic and boundary integers c, converted at the working precision P
// as in multimodular GEMM: the BLAS CRT must be within its bound
//   |v - c| < 2 num_primes 2^21 2^(24 t0) + 2^-P |c|
// (blas_residues.hxx). It is not always bitwise equal to the exact CRT:
// c = 2^K is the standard counterexample (the error borrows across the kept
// limbs), and c = 1 lies within the CRT error of zero.
TEST_CASE("blas_residues: exact dyadic CRT cases")
{
  if(El::mpi::Rank() != 0)
    return;
  // Precision of the mpf results; the Context is that of the multimodular
  // Cholesky/Trsm at this precision (bits = P + 64), independent of
  // El::gmp::Precision()
  const long P = GENERATE(long(768), long(1024));
  const long bits = P + 64;
  const slong k = 1;
  DYNAMIC_SECTION("P=" << P)
  {
    const auto &context = Blas_Residues::get_context(bits, bits, k);
    const size_t num_primes = context.num_primes;
    CAPTURE(num_primes);
    CAPTURE(context.t0);
    REQUIRE(context.t0 > 0);

    // c = s 2^K + d, for 2^K at the top of the product range, at a limb
    // boundary, in the middle and just above the CRT error; and small c
    const long K_max = 2 * bits - 2;
    std::vector<long> exponents{K_max, 2 * bits - 64 - 1, 2 * bits - 128,
                                bits,
                                Blas_Residues::chunk_bits * context.t0 + 64};
    std::vector<fmpz> c;
    auto add = [&c](const fmpz_t x) {
      c.emplace_back();
      fmpz_init_set(&c.back(), x);
    };
    fmpz_t x;
    fmpz_init(x);
    for(const long K : exponents)
      for(const int sign : {1, -1})
        for(const int d : {0, 1, -1})
          {
            fmpz_one(x);
            fmpz_mul_2exp(x, x, K);
            if(sign < 0)
              fmpz_neg(x, x);
            if(d > 0)
              fmpz_add_ui(x, x, d);
            else if(d < 0)
              fmpz_sub_ui(x, x, -d);
            add(x);
          }
    for(const long small : {0L, 1L, -1L, 3L})
      {
        fmpz_set_si(x, small);
        add(x);
      }
    fmpz_clear(x);
    const size_t n = c.size();

    std::vector<int32_t> y(n * num_primes), y_flint(n * num_primes);
    std::vector<double> quotients(n, 0.0);
    for(size_t p = 0; p < num_primes; ++p)
      {
        const auto prime = context.comb.primes[p];
        for(size_t e = 0; e < n; ++e)
          {
            const ulong r = fmpz_fdiv_ui(&c[e], prime);
            y_flint[e + n * p] = r;
            const double r_u = double((r * ulong(context.u[p])) % prime);
            y[e + n * p] = int32_t(r_u);
            quotients[e] += r_u * context.inverse_primes[p];
          }
      }

    Blas_Residues::Crt_Workspace work(context);
    Blas_Residues::blas_crt_high_chunks(context, y.data(), n, n, work);
    mpf_t blas_value, flint_value, exact, diff, bound, absolute_bound;
    mpf_init2(blas_value, P);
    mpf_init2(flint_value, P);
    const mp_bitcnt_t exact_precision = 2 * bits + 128;
    for(auto *v : {&exact, &diff, &bound, &absolute_bound})
      mpf_init2(*v, exact_precision);
    // 2 num_primes 2^21 2^(24 t0)
    mpf_set_ui(absolute_bound, 2 * num_primes);
    mpf_mul_2exp(absolute_bound, absolute_bound,
                 21 + Blas_Residues::chunk_bits * context.t0);
    for(size_t e = 0; e < n; ++e)
      {
        CAPTURE(e);
        CAPTURE(fmpz_bits(&c[e]));
        mpf_set_ui(blas_value, 0);
        mpf_set_ui(flint_value, 0);
        Blas_Residues::blas_crt_to_mpf(context, work.S.data() + e, n,
                                       quotients[e], 0, work, blas_value);
        Blas_Residues::flint_crt_to_mpf(context, y_flint.data() + e, n, 0,
                                        work, flint_value);
        mpz_t z;
        mpz_init(z);
        fmpz_get_mpz(z, &c[e]);
        mpf_set_z(exact, z);
        mpz_clear(z);

        // 2^-P |c|
        mpf_abs(bound, exact);
        mpf_div_2exp(bound, bound, P);
        // FLINT: c truncated to P bits
        mpf_sub(diff, flint_value, exact);
        mpf_abs(diff, diff);
        REQUIRE(mpf_cmp(diff, bound) <= 0);
        // BLAS: within the CRT bound
        mpf_add(bound, bound, absolute_bound);
        mpf_sub(diff, blas_value, exact);
        mpf_abs(diff, diff);
        CAPTURE(mpf_get_d(diff));
        REQUIRE(mpf_cmp(diff, bound) < 0);
      }
    for(auto *v :
        {&blas_value, &flint_value, &exact, &diff, &bound, &absolute_bound})
      mpf_clear(*v);
    for(auto &v : c)
      fmpz_clear(&v);
  }
}
