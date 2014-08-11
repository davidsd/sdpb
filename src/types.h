#ifndef SDP_BOOTSTRAP_TYPES_H_
#define SDP_BOOTSTRAP_TYPES_H_

#include <mblas.h>
#include <mlapack.h>

typedef mpackint Integer;
typedef mpf_class Real;

double realToDouble(Real r) {
  return mpf_get_d(r.get_mpf_t());
}

#endif  // SDP_BOOTSTRAP_TYPES_H_
