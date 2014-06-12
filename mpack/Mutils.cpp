/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Mutils.cpp,v 1.7 2009/09/16 08:32:46 nakatamaho Exp $ 
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

#include <mblas_gmp.h>
#include <mlapack_gmp.h>
#include <math.h>
#define ___MPACK_BUILD_WITH_GMP___

mpf_class
mpf_approx_log2(mpf_class x)
{
#if defined ___MPACK_BUILD_WITH_GMP___
    double d;
    double ln2_app;
    signed long int exp;

    d = mpf_get_d_2exp(&exp, x.get_mpf_t());
    ln2_app = (double)exp + log10(d) / log10(2);
    return ln2_app;
#endif
#if defined ___MPACK_BUILD_WITH_QD___
  return log10(x) / (qd_real::_log2/qd_real::_log10);
#endif
#if defined ___MPACK_BUILD_WITH_DD___
  return log10(x) / (dd_real::_log2/dd_real::_log10);
#endif
}

mpf_class
mpf_approx_log(mpf_class x)
{
#if defined ___MPACK_BUILD_WITH_GMP___
    double d;
    double ln_app;
    signed long int exp;

    d = mpf_get_d_2exp(&exp, x.get_mpf_t());
    ln_app = (double)exp * log (2.0) + log(d);
    return ln_app;
#endif
#if defined ___MPACK_BUILD_WITH_QD___
    return log(x);
#endif
#if defined ___MPACK_BUILD_WITH_DD___
    return log(x);
#endif
}

mpf_class
mpf_approx_log10(mpf_class x)
{
#if defined ___MPACK_BUILD_WITH_GMP___
    double d;
    double ln10_app;
    signed long int exp;

    d = mpf_get_d_2exp(&exp, x.get_mpf_t());
    ln10_app = (double)exp * log10(2.0) + log10(d);
    return ln10_app;
#endif
#if defined ___MPACK_BUILD_WITH_QD___
    return log10(x);
#endif
#if defined ___MPACK_BUILD_WITH_DD___
    return log10(x);
#endif
}

mpf_class
mpf_approx_pow(mpf_class x, mpf_class y)
{
#if defined ___MPACK_BUILD_WITH_GMP___
    mpf_class mtemp1, mtemp2;
    mtemp1 = y * mpf_approx_log(x);
    mtemp2 = mpf_approx_exp(mtemp1);
    return mtemp2;
#endif
#if defined ___MPACK_BUILD_WITH_QD___
    return pow(x, y);
#endif
#if defined ___MPACK_BUILD_WITH_DD___
    return pow(x, y);
#endif
}

mpf_class
mpf_approx_cos(mpf_class x)
{
#if defined ___MPACK_BUILD_WITH_GMP___
    mpf_class mtemp1;
    mtemp1 = cos(x.get_d());
    return mtemp1;
#endif
#if defined ___MPACK_BUILD_WITH_QD___
    return cos(x);
#endif
#if defined ___MPACK_BUILD_WITH_DD___
    return cos(x);
#endif
}

mpf_class
mpf_approx_sin(mpf_class x)
{
#if defined ___MPACK_BUILD_WITH_GMP___
    mpf_class mtemp1;
    mtemp1 = sin(x.get_d());
    return mtemp1;
#endif
#if defined ___MPACK_BUILD_WITH_QD___
    return sin(x);
#endif
#if defined ___MPACK_BUILD_WITH_DD___
    return sin(x);
#endif
}

mpf_class
mpf_approx_exp(mpf_class x)
{
#if defined ___MPACK_BUILD_WITH_GMP___
    mpf_class mtemp1;
    mtemp1 = exp(x.get_d());
    return mtemp1;
#endif
#if defined ___MPACK_BUILD_WITH_QD___
    return exp(x);
#endif
#if defined ___MPACK_BUILD_WITH_DD___
    return exp(x);
#endif
}

mpf_class
mpf_approx_pi()
{
#if defined ___MPACK_BUILD_WITH_GMP___
    mpf_class mtemp1;
    mtemp1 = M_PI;
    return mtemp1;
#endif
#if defined ___MPACK_BUILD_WITH_QD___
    return qd_real::_pi;
#endif
#if defined ___MPACK_BUILD_WITH_DD___
    return dd_real::_pi;
#endif
}
