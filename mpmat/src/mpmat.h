//
// Created by Petr Kravchuk on 8/10/17.
//
// *****************
// *    Warning    *
// *****************
//
// We use a lot of GMP internals, and thus may be incompatible with future versions
//    Last tested version : GMP 6.1.2
//


#ifndef MPMAT_MPMAT_H
#define MPMAT_MPMAT_H

#include <gmpxx.h>

// The floating point type used to emulate integer arithmetic
typedef double mpmat_double;

// The effective number of bits in a mpmat_double mantissa
// The actual number of bits (explicitly stored) for mantissa can be lower
// We assume that mpmat_double can represent MPMAT_DOUBLE_MANT_IMPLICIT-bit integer value without loss of information
// These bit counts do not include sign bit
#define MPMAT_DOUBLE_MANT_IMPLICIT 53

// Converts a single mpf_class into an array of mpmat_double
//
// Arguments:
// (in)  source    : the mpf_class to convert
// (out) dest      : pointer to the start of the destination; must be allocated
// (in)  size      : size of the space allocated at dest, in units of mpmat_double
// (in)  mpmat_limb: the limb size to use inside mpmat_double; must be <= MPMAT_DOUBLE_MANT_IMPLICIT
// (in)  exp       : the exponent for the resulting array
//
// * dest[0] will contain the most significant limb of the number (this is reversed relative to GMP)
// * exp is the exponent, in bits, of the resulting number
// * dest is filled completely, either filling the least significant bits with 0 or discarding those of GMP,
//   depending on the relation between size and _mp_size
// * mpmat_limb should be no greater than the number of bits in mp_limb_t
// * we assume that mpmat_double 0 is represented by a zero bitstring (i.e. IEEE754)
//
void mpmatConvertGMPToDouble(const mpf_class source,
                             mpmat_double * dest,
                             const int size,
                             const int mpmat_limb,
                             const int exp);
// TODO:                     const int stride=1); // in order to distribute into different matrices
// TODO: If these conversion steps become critical, we may check whether we miss a lot of cache (we probably do, but who
// TODO:   knows how clever the prefetching is). In a non-striding scenario we need a cache-friendly transposition code.
// TODO: A better idea is probably to convert without transposition and then transpose using an efficient BLAS func

// Converts an array of mpmat_double into mpf_class
//
// Arguments:
// (out) dest      : the destination mpf_variable
// (in)  source    : pointer the start of the mpmat_double array
// (in)  size      : length of source
// (in)  mpmat_limb: the original limb size of the source
// (in)  exp       : the exponent of the array, in bits
//
// * mpmat_double in source can contain more bits then mpmat_limb, up to MPMAT_DOUBLE_MANT_IMPLICIT, they are extracted
//   by static_cast<mp_limb_t>; mpmat_limb essentially gives the relative exponent of different mpmat_double's
// * This function fills all limbs available in dest, padding by 0 or discarding the least significant limbs
// * We assume that MPMAT_DOUBLE_MANT_IMPLICIT < 8*sizeof(mp_limb_t) and that mpmat_limb > 0 (strict inequalities)
// * If some mpmat_double are negative, it may happen that some of the leading limbs of the result, which would be
//   non-zero if all mpmat_double were positive and occupied MPMAT_DOUBLE_MANT_IMPLICIT bits, vanish. In this case the
//   result is correctly shifted in order to produce a valid mpf_t, but the lower limbs which were discarded during the
//   calculation, are not restored. This is not a problem if the precision of dest is enough to fit all limbs of source,
//   and no information is ever discarded.
//   More generally, because to various exp-related paddings, one limb of precision on dest can be lost currently.
//   TODO: do we need to change this behavior?
//
void mpmatConvertDoubleToGMP(mpf_class & dest,
                             const mpmat_double * source,
                             const int size,
                             const int mpmat_limb,
                             const int exp);
// TODO:                     const int stride=1); // in order to read from different matrices

// Performs multilplication of GMP numbers a and b by converting to mpmat_double arrays, and stores result in dest.
// Doesn't loose precision except as described in mpmatConvertDoubleToGMP
// Test function before going to matrices
//
// Arguments:
// (out) dest : initialized mpf_class to store the result
// (in)  a    : first operand
// (in)  b    : second operand
//
void mpmatMultiplyGMPBaseCase(mpf_class & dest,
                              const mpf_class  a,
                              const mpf_class  b);


#endif //MPMAT_MPMAT_H
