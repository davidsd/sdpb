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
// * !!!!!! Currently contains a bug, see code
void mpmatConvertDoubleToGMP(mpf_class & dest,
                             const mpmat_double * source,
                             const int size,
                             const int mpmat_limb,
                             const int exp);

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

// Converts GMP vector to array of mpmat_double vectors; all operations are performed in-place over dest
//
// Arguments:
// (in)  source     : the pointer to the first element in the GMP vector
// (in)  source_len : the length of the GMP vector
// (out) dest       : the pointer to the allocated array of mpmat_doubles
// (in)  size       : number of mpmat_doubles to use per GMP entry
// (in)  mpmat_limb : the bit limb size to use in mpmat_double array
// (out) exp        : the exponent used in the output
//
// * The same exponent exp is used for all entries of the vector
// * This function automatically detects the maximal exponent in source and sets exp accordingly
// * This function internally uses mpmatConvertGMPToDouble
//
void mpmatConvertGMPToDoubleVector(const mpf_class * source,
                                   const int source_len,
                                   mpmat_double * dest,
                                   const int mpmat_size,
                                   const int mpmat_limb,
                                   int & exp);

// Converts an array of mpmat_double vectors into a GMP vector; in-place transforms the mpmat_double array
//
// Arguments:
// (out)    dest       : the destination GMP vector; must be allocated & mpf's intialized
// (in)     dest_len   : the length of the vector
// (in/out) source     : the array of source vectors
// (in)     mpmat_size : size of mpmat_double arrays
// (in)     mpmat_limb : the original limb size of mpmat_double
// (in)     exp        : the exponent to use for interpretation of source
//
// * Note that source is transposed in-place after the call
// * This function internally uses mpmatConvertDoubleToGMP
//
void mpmatConvertDoubleToGMPVector(mpf_class * dest,
                                   const int dest_len,
                                   mpmat_double * source,
                                   const int mpmat_size,
                                   const int mpmat_limb,
                                   int exp);
#endif //MPMAT_MPMAT_H
