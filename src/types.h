#ifndef SDP_BOOTSTRAP_TYPES_H_
#define SDP_BOOTSTRAP_TYPES_H_

#include <mblas.h>
#include <mlapack.h>

#if defined ___MPACK_BUILD_WITH_GMP___
typedef mpackint Integer;
typedef mpf_class Real;

inline void setDefaultPrecision(int prec) {
  mpf_set_default_prec(prec);
}

inline int getDefaultPrecision() {
  return mpf_get_default_prec();
}

inline void setPrecision(Real &r, int prec) {
  r.set_prec(prec);
}

inline int getPrecision(const Real &r) {
  return r.get_prec();
}
#endif

#if defined ___MPACK_BUILD_WITH_MPFR___
typedef mpackint Integer;
typedef mpreal Real;

inline void setDefaultPrecision(int prec) {
  mpreal::set_default_prec(prec);
}

inline int getDefaultPrecision() {
  return mpreal::get_default_prec();
}

inline void setPrecision(Real &r, int prec) {
  r.set_prec(prec);
}

inline int getPrecision(const Real &r) {
  return r.get_prec();
}
#endif

#endif  // SDP_BOOTSTRAP_TYPES_H_
