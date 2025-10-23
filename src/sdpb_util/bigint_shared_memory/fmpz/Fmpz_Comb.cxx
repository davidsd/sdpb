#include "Fmpz_Comb.hxx"

#include "sdpb_util/assert.hxx"

#include <flint/ulong_extras.h>

#include <El.hpp>

// TODO explain all parameters!
namespace
{
#define MAX_BLAS_DP_INT (UWORD(1) << 53)

  mp_limb_t calculate_output_bits(mp_limb_t Abits, mp_limb_t Bbits, int sign,
                                  mp_limb_signed_t k)
  {
    return Abits + Bbits + FLINT_BIT_COUNT(k) + sign;
  }

  // adapted from FLINT, fmpz_mat/mul_blas.c
  // TODO replace C-style allocation with std::vector::resize()
  mp_limb_t *
  _calculate_primes(slong *num_primes_, flint_bitcnt_t bits, slong k)
  {
    slong num_primes, primes_alloc;
    mp_limb_t *primes;
    mp_limb_t p;
    fmpz_t prod;

    p = 2 + 2 * n_sqrt((MAX_BLAS_DP_INT - 1) / (ulong)k);
    if(bits > 200)
      {
        /* if mod is the bottleneck, ensure p1*p2*p3 < 2^62 */
        p = FLINT_MIN(p, UWORD(1664544));
      }

    primes_alloc = 1 + bits / FLINT_BIT_COUNT(p);
    primes = FLINT_ARRAY_ALLOC(primes_alloc, mp_limb_t);
    num_primes = 0;

    fmpz_init_set_ui(prod, 1);

    do
      {
        do
          {
            if(p < 1000)
              {
                fmpz_clear(prod);
                flint_free(primes);
                *num_primes_ = 0;
                return NULL;
              }
            p--;
        } while(!n_is_prime(p));

        if(num_primes + 1 > primes_alloc)
          {
            primes_alloc = FLINT_MAX(num_primes + 1, primes_alloc * 5 / 4);
            primes = FLINT_ARRAY_REALLOC(primes, primes_alloc, mp_limb_t);
          }

        primes[num_primes] = p;
        num_primes++;

        fmpz_mul_ui(prod, prod, p);

    } while(fmpz_bits(prod) <= bits);

    fmpz_clear(prod);

    *num_primes_ = num_primes;
    return primes;
  }

  std::vector<mp_limb_t> calculate_primes(flint_bitcnt_t bits, slong k)
  {
    slong n;
    mp_limb_t *primes = _calculate_primes(&n, bits, k);
    // TODO test this case (happens with low precision) and add workaround
    ASSERT(primes != NULL, "Failed to calculate primes for bits=", bits,
           "k=", k);

    auto result = std::vector<mp_limb_t>(primes, primes + n);
    flint_free(primes);
    return result;
  }
}

Fmpz_Comb::Fmpz_Comb(mp_limb_t bits, mp_limb_signed_t k)
    : primes(calculate_primes(bits, k)),
      num_primes(primes.size()),
      mods(num_primes),
      shifts(num_primes)
{
  fmpz_comb_init(comb, primes.data(), num_primes);
  fmpz_comb_temp_init(comb_temp, comb);
  for(size_t i = 0; i < num_primes; ++i)
    {
      auto &mod = mods.at(i);
      nmod_init(&mod, primes.at(i));
      shifts.at(i) = ((2 * MAX_BLAS_DP_INT) / mod.n) * mod.n;
    }
}

Fmpz_Comb::Fmpz_Comb(mp_limb_t Abits, mp_limb_t Bbits, int sign,
                     mp_limb_signed_t k)
    : Fmpz_Comb(calculate_output_bits(Abits, Bbits, sign, k), k)
{}

Fmpz_Comb::~Fmpz_Comb()
{
  fmpz_comb_temp_clear(comb_temp);
  fmpz_comb_clear(comb);
}
