/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rsytd2.cpp,v 1.6 2009/09/22 21:28:58 nakatamaho Exp $ 
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
Rsytd2(const char *uplo, mpackint n, mpf_class * A, mpackint lda, mpf_class * d,
    mpf_class * e, mpf_class * tau, mpackint *info)
{

    mpf_class One = 1.0, Zero = 0.0, Half = 0.5;
    mpf_class taui, alpha;
    mpackint upper;
    mpackint i;

    *info = 0;
    upper = Mlsame_gmp(uplo, "U");
    if (!upper && !Mlsame_gmp(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((mpackint)1, n)) {
	*info = -4;
    }
    if (*info != 0) {
	Mxerbla_gmp("Rsytd2", -(*info));
	return;
    }
//Quick return if possible
    if (n <= 0)
	return;
    if (upper) {
//Reduce the upper triangle of A
	for (i = n - 1; i >= 1; i--) {
//Generate elementary reflector H(i) = I - tau * v * v'
//to annihilate A(1:i-1,i+1)
	    Rlarfg(i, &A[(i - 1) + i * lda], &A[0 + i * lda], 1, &taui);
	    e[i - 1] = A[(i - 1) + i * lda];
	    if (taui != Zero) {
//Apply H(i) from both sides to A(1:i,1:i)
		A[(i - 1) + i * lda] = One;
//Compute  x := tau * A * v  storing x in TAU(1:i)
		Rsymv(uplo, i, taui, A, lda, &A[0 + i * lda], 1, Zero, tau, 1);
//Compute  w := x - 1/2 * tau * (x'*v) * v
		alpha = -Half * taui * Rdot(i, tau, 1, &A[0 + i * lda], 1);
		Raxpy(i, alpha, &A[0 + i * lda], 1, tau, 1);
//Apply the transformation as a rank-2 update
//A := A - v * w' - w * v'
		Rsyr2(uplo, i, -One, &A[0 + i * lda], 1, tau, 1, A, lda);
		A[(i - 1) + i * lda] = e[i - 1];
	    }
	    d[i] = A[i + i * lda];
	    tau[i - 1] = taui;
	}
	d[0] = A[0];
    } else {
//Reduce the lower triangle of A
	for (i = 1; i <= n - 1; i++) {
//Generate elementary reflector H(i) = I - tau * v * v'
//to annihilate A(i+2:n,i)
	    Rlarfg(n - i, &A[i + (i - 1) * lda], &A[min(i + 2,
			n) - 1 + (i - 1) * lda], 1, &taui);
	    e[i - 1] = A[i + (i - 1) * lda];
	    if (taui != Zero) {
//Apply H(i) from both sides to A(i+1:n,i+1:n)
		A[i + (i - 1) * lda] = One;
//Compute  x := tau * A * v  storing y in TAU(i:n-1)
		Rsymv(uplo, n - i, taui, &A[i + i * lda],
		    lda, &A[i + (i - 1) * lda], 1, Zero, &tau[i - 1], 1);
//Compute  w := x - 1/2 * tau * (x'*v) * v
		alpha =
		    -Half * taui * Rdot(n - i, &tau[i - 1], 1,
		    &A[i + (i - 1) * lda], 1);
		Raxpy(n - i, alpha, &A[i + (i - 1) * lda], 1, &tau[i - 1], 1);
//Apply the transformation as a rank-2 update:
//A := A - v * w' - w * v'
		Rsyr2(uplo, n - i, -One, &A[i + (i - 1) * lda], 1, &tau[i - 1],
		    1, &A[i + i * lda], lda);
		A[i + (i - 1) * lda] = e[i - 1];
	    }
	    d[i - 1] = A[(i - 1) + (i - 1) * lda];
	    tau[i - 1] = taui;
	}
	d[n - 1] = A[(n - 1) + (n - 1) * lda];
    }
    return;
}
