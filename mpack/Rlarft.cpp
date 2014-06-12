/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rlarft.cpp,v 1.2 2009/09/12 07:59:10 nakatamaho Exp $ 
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

#include <mblas_gmp.h>
#include <mlapack_gmp.h>

void
Rlarft(const char *direct, const char *storev, mpackint n, mpackint k, mpf_class * v,
    mpackint ldv, mpf_class * tau, mpf_class * t, mpackint ldt)
{
    mpf_class Zero = 0.0, One = 1.0;
    mpf_class vii;
    mpackint i, j;

    //Quick return if possible
    if (n == 0)
	return;

    if (Mlsame_gmp(direct, "F")) {
	for (i = 1; i <= k; i++) {
	    if (tau[i - 1] == Zero) {
		//H(i)  =  I
		for (j = 1; j <= i; j++) {
		    t[(j - 1) + (i - 1) * ldt] = Zero;
		}
	    } else {
		//general case
		vii = v[(i - 1) + (i - 1) * ldv];
		v[(i - 1) + (i - 1) * ldv] = One;
		if (Mlsame_gmp(storev, "C")) {
		    // T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
		    Rgemv("Transpose", n - i + 1, i - 1, -tau[i - 1],
			&v[(i - 1) + 0 * ldv], ldv,
			&v[(i - 1) + (i - 1) * ldv], 1, Zero,
			&t[0 + (i - 1) * ldt], 1);
		} else {
		    //T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
		    Rgemv("No transpose", i - 1, n - i + 1, -tau[i - 1],
			&v[0 + (i - 1) * ldv], ldv,
			&v[(i - 1) + (i - 1) * ldv], ldv, Zero,
			&t[0 + (i - 1) * ldt], 1);
		}
		v[(i - 1) + (i - 1) * ldv] = vii;
		//T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
		Rtrmv("Upper", "No transpose", "Non-unit", i - 1, t, ldt,
		    &t[0 + (i - 1) * ldt], 1);
		t[(i - 1) + (i - 1) * ldt] = tau[i - 1];
	    }
	}
    } else {
	for (i = k; i >= 1; i--) {
	    if (tau[i - 1] == Zero) {
		//H(i)  =  I
		for (j = i; j < k; j++) {
		    t[(j - 1) + (i - 1) * ldt] = Zero;
		}
	    } else {
		//general case
		if (i < k) {
		    if (Mlsame_gmp(storev, "C")) {
			vii = v[(n - k + i - 1) + (i - 1) * ldv];
			v[(n - k + i - 1) + (i - 1) * ldv] = One;
			//T(i+1:k,i) := - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
			Rgemv("Transpose", n - k + i, k - i, -tau[i - 1],
			    &v[0 + i * ldv], ldv, &v[0 + (i - 1) * ldv], 1,
			    Zero, &t[i + (i - 1) * ldt], 1);
			v[(n - k + i - 1) + (i - 1) * ldv] = vii;
		    } else {
			vii = v[(i - 1) + (n - k + i - 1) * ldv];
			v[(i - 1) + (n - k + i - 1) * ldv] = One;
			//T(i+1:k,i) := - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
			Rgemv("No transpose", k - i, n - k + i, -tau[i - 1],
			    &v[i + 0 * ldv], ldv, &v[(i - 1) + 0 * ldv], ldv,
			    Zero, &t[i + (i - 1) * ldt], 1);
			v[(i - 1) + (n - k + i - 1) * ldv] = vii;
		    }
		    //T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
		    Rtrmv("Lower", "No transpose", "Non-unit", k - i,
			&t[i + i * ldt], ldt, &t[i + (i - 1) * ldt], 1);
		}
		t[(i - 1) + (i - 1) * ldt] = tau[i - 1];
	    }
	}
    }
    return;
}
