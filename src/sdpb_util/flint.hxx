#pragma once

// In 3.2.0, flint.h does not include gmp.h
// see https://github.com/flintlib/flint/pull/1969
// It should be included before fmpz.h to enable fmpz_get_mpf()
#include <gmp.h>
#include <flint/fmpz.h>
// In FLINT 2.8.5, part of nmod_vec.h was extracted to nmod.h
// See:
// https://github.com/davidsd/sdpb/issues/234
// https://github.com/flintlib/flint/pull/1041
#if __FLINT_RELEASE >= 20805
#include <flint/nmod.h>
#else
#include <flint/nmod_vec.h>
#endif
