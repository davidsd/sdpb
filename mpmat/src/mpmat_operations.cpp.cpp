//
// Created by Petr Kravchuk on 8/14/17.
//

#include "mpmat.h"
#include <gmpxx.h>
#include <math.h>
#include <cassert>
#include <iostream>

template <typename T>
inline T ceil_div(T a, T b) {
    return a / b + ( (a % b != 0) & (a > 0) );
}

void mpmatMultiplyGMPBaseCase(mpf_class & dest,
                              const mpf_class  a,
                              const mpf_class  b) {

    int mpmat_limb = MPMAT_DOUBLE_MANT_IMPLICIT / 2;
    int mpmat_size_a = ceil_div( a.get_mpf_t()->_mp_size * mp_bits_per_limb, mpmat_limb );
    int mpmat_size_b = ceil_div( b.get_mpf_t()->_mp_size * mp_bits_per_limb, mpmat_limb );

    while ( 2 * mpmat_limb + ceil(log2(fmax(mpmat_size_a, mpmat_size_b))) > MPMAT_DOUBLE_MANT_IMPLICIT ) {
        mpmat_limb = ( MPMAT_DOUBLE_MANT_IMPLICIT - ceil(log2(fmax(mpmat_size_a, mpmat_size_b))) ) / 2;
        mpmat_size_a = ceil_div( a.get_mpf_t()->_mp_size * mp_bits_per_limb, mpmat_limb );
        mpmat_size_b = ceil_div( b.get_mpf_t()->_mp_size * mp_bits_per_limb, mpmat_limb );
    }

    int mpmat_size_c = fmax(mpmat_size_a, mpmat_size_b);

    assert(mpmat_limb > 0);

    mpmat_double * a_double = new mpmat_double[mpmat_size_a];
    mpmat_double * b_double = new mpmat_double[mpmat_size_b];
    mpmat_double * c_double = new mpmat_double[mpmat_size_c];

    int a_mp_exp = a.get_mpf_t() -> _mp_exp * mp_bits_per_limb;
    int b_mp_exp = b.get_mpf_t() -> _mp_exp * mp_bits_per_limb;

    mpmatConvertGMPToDouble(a,a_double,mpmat_size_a, mpmat_limb, a_mp_exp);
    mpmatConvertGMPToDouble(b,b_double,mpmat_size_b, mpmat_limb, b_mp_exp);

    for(int i = 0; i< mpmat_size_c; i++) {
        c_double[i]=0;
        for(int k = 0; k<=i; k++) {
            c_double[i] += a_double[i-k] * b_double[k];
        }
    }

    mpmatConvertDoubleToGMP(dest, c_double, mpmat_size_c,mpmat_limb, a_mp_exp+b_mp_exp-mpmat_limb);

    delete [] a_double;
    delete [] b_double;
    delete [] c_double;
}
