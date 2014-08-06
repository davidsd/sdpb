#ifndef SDP_BOOTSTRAP_TYPES_H_
#define SDP_BOOTSTRAP_TYPES_H_

#include <mblas_gmp.h>
#include <mlapack_gmp.h>

typedef mpz_class Integer;
typedef mpf_class Real;

double realToDouble(Real r) {
  return mpf_get_d(r.get_mpf_t());
}

#endif  // SDP_BOOTSTRAP_TYPES_H_
