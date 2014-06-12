/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rlaev2.cpp,v 1.4 2009/09/26 02:21:32 nakatamaho Exp $ 
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
/*
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer. 
  
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.
  
- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.
  
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
*/

// http://www.netlib.org/lapack/double/dlaev2.f

#include <mblas_gmp.h>
#include <mlapack_gmp.h>
#include <stdio.h> //for printf. shall be removed

void
Rlaev2(mpf_class a, mpf_class b, mpf_class c, mpf_class * rt1, mpf_class * rt2,
    mpf_class * cs1, mpf_class * sn1)
{
    mpf_class ab, acmn, acmx, acs, adf;
    mpf_class cs, ct, df, rt, sm, tb, tn;
    mpf_class zero, one, two, half;
    mpackint sgn1, sgn2;

    zero = 0.0;
    one = 1.0;
    two = 2.0;
    half = 0.5;

    sm = a + c;
    df = a - c;
    adf = abs(df);
    tb = b + b;
    ab = abs(tb);

    if (abs(a) > abs(c)) {
	acmx = a;
	acmn = c;
    } else {
	acmx = c;
	acmn = a;
    }
    if (adf > ab) {
	rt = adf * sqrt(one + (ab / adf) * (ab / adf));
    } else if (adf < ab) {
	rt = ab * sqrt(one + (adf / ab) * (adf / ab));
    } else {
//Includes case AB=ADF=0
	rt = ab * sqrt(two);
    }
    if (sm < zero) {
	*rt1 = half * (sm - rt);
	sgn1 = -1;
//Order of execution important.
//To get fully accurate smaller eigenvalue,
//next line needs to be executed in higher precision.
	*rt2 = (acmx / (*rt1)) * acmn - (b / (*rt1)) * b;
    } else if (sm > zero) {
	*rt1 = half * (sm + rt);
	sgn1 = 1;
//Order of execution important.
//To get fully accurate smaller eigenvalue,
//next line needs to be executed in higher precision.
	*rt2 = (acmx / (*rt1)) * acmn - (b / (*rt1)) * b;
    } else {
//Includes case RT1 = RT2 = 0
	*rt1 = half * rt;
	*rt2 = -1.0 * half * rt;
	sgn1 = 1;
    }
//Compute the eigenvector
    if (df >= zero) {
	cs = df + rt;
	sgn2 = 1;
    } else {
	cs = df - rt;
	sgn2 = -1;
    }
    acs = abs(cs);
    if (acs > ab) {
	ct = -tb / cs;
	*sn1 = one / sqrt(one + ct * ct);
	*cs1 = ct * (*sn1);
    } else {
	if (ab == zero) {
	    *cs1 = one;
	    *sn1 = zero;
	} else {
	    printf("#Rlaev2 Checkpoint 13 Not checked\n");
            exit(1);
	    tn = -cs / tb;
	    *cs1 = one / sqrt(one + tn * tn);
	    *sn1 = tn * (*cs1);
	}
    }
    if (sgn1 == sgn2) {
	tn = *cs1;
	*cs1 = -(*sn1);
	*sn1 = tn;
    }
    return;
}
