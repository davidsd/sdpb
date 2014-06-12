/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2009 by Nakata, Maho
 *
 * MPACK - multiple precision arithmetic library
 *
 * This file is part of MPACK.
 *
 * MPACK is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License version 3
 * only, as published by the Free Software Foundation.
 *
 * MPACK is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License version 3 for more details
 * (a copy is included in the LICENSE file that accompanied this code).
 *
 * You should have received a copy of the GNU Lesser General Public License
 * version 3 along with MPACK.  If not, see
 * <http://www.gnu.org/licenses/lgpl.html>
 * for a copy of the LGPLv3 License.
 *
 ************************************************************************/

#ifndef _MUTILS_GMP_H_
#define _MUTILS_GMP_H_

using std::max;
using std::min;

mpf_class Msign(mpf_class a, mpf_class b);
double cast2double(mpf_class a);
int M2int(mpf_class a);
void mpf_pow(mpf_t ans, mpf_t x, mpf_t y);
mpf_class mpf_approx_log(mpf_class x);
mpf_class mpf_approx_log2(mpf_class x);
mpf_class mpf_approx_log10(mpf_class x);
mpf_class mpf_approx_pow(mpf_class x, mpf_class y);
mpf_class mpf_approx_cos(mpf_class x);
mpf_class mpf_approx_sin(mpf_class x);
mpf_class mpf_approx_exp(mpf_class x);
mpf_class mpf_approx_pi();

//implementation of sign transfer function.
inline mpf_class
Msign(mpf_class a, mpf_class b)
{
    mpf_class mtmp;
    mpf_abs(mtmp.get_mpf_t(), a.get_mpf_t());
    if (b != 0.0) {
	mtmp = mpf_sgn(b.get_mpf_t()) * mtmp;
    }
    return mtmp;
}

inline double
cast2double(mpf_class a)
{
    return a.get_d();
}

inline int
M2int(mpf_class a)
{
    int i;
    mpf_t tmp;
    a = a + 0.5;
    mpf_floor(tmp, a.get_mpf_t());
    i = (int)mpf_get_si(tmp);
    return i;
}

#endif
