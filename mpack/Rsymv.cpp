/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rsymv.cpp,v 1.4 2009/09/24 07:25:57 nakatamaho Exp $ 
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

/*
Based on http://www.netlib.org/blas/dsymv.f
Rsymv performs the matrix-vector  operation
 y := alpha*A*x + beta*y,
where alpha and beta are scalars, x and y are n element vectors and
 A is an n by n symmetric matrix.
*/

#include <mblas_gmp.h>

void
Rsymv(const char *uplo, mpackint n, mpf_class alpha, mpf_class * A,
    mpackint lda, mpf_class * x, mpackint incx, mpf_class beta, mpf_class * y,
    mpackint incy)
{
    mpackint ix, iy, jx, jy, kx, ky;

    mpf_class Zero = 0.0, One = 1.0;

    mpf_class temp1, temp2;

    //test the input parameters.
    mpackint info = 0;

    if (!Mlsame_gmp(uplo, "U") && !Mlsame_gmp(uplo, "L"))
	info = 1;
    else if (n < 0)
	info = 2;
    else if (lda < max((mpackint) 1, n))
	info = 5;
    else if (incx == 0)
	info = 7;
    else if (incy == 0)
	info = 10;

    if (info != 0) {
	Mxerbla_gmp("Rsymv ", info);
	return;
    }
    //quick return if possible.
    if ((n == 0) || ((alpha == Zero) && (beta == One)))
	return;

    //set up the start points in  x  and  y.
    if (incx > 0)
	kx = 0;
    else
	kx = -(n - 1) * incx;
    if (incy > 0)
	ky = 0;
    else
	ky = -(n - 1) * incy;

    //start the operations. in this version the elements of a are
    //accessed sequentially with one pass through the triangular part
    //of A.

    //first form  y := beta*y.
    if (beta != One) {
	iy = ky;
	if (beta == Zero) {
	    for (mpackint i = 0; i < n; i++) {
		y[iy] = Zero;
		iy = iy + incy;
	    }
	} else {
	    for (mpackint i = 0; i < n; i++) {
		y[iy] = beta * y[iy];
		iy = iy + incy;
	    }
	}
    }
    if (alpha == Zero)
	return;

    if (Mlsame_gmp(uplo, "U")) {
	//form  y  when a is stored in upper triangle.
	jx = kx;
	jy = ky;
	for (mpackint j = 0; j < n; j++) {
	    temp1 = alpha * x[jx];
	    temp2 = Zero;
	    ix = kx;
	    iy = ky;
	    for (mpackint i = 0; i < j; i++) {
		y[iy] = y[iy] + temp1 * A[i + j * lda];
		temp2 = temp2 + A[i + j * lda] * x[ix];
		ix = ix + incx;
		iy = iy + incy;
	    }
	    y[jy] = y[jy] + temp1 * A[j + j * lda] + alpha * temp2;
	    jx = jx + incx;
	    jy = jy + incy;
	}
    } else {
	//form  y  when a is stored in lower triangle.
	jx = kx;
	jy = ky;
	for (mpackint j = 0; j < n; j++) {
	    temp1 = alpha * x[jx];
	    temp2 = Zero;
	    y[jy] = y[jy] + temp1 * A[j + j * lda];
	    ix = jx;
	    iy = jy;
	    for (mpackint i = j + 1; i < n; i++) {
		ix = ix + incx;
		iy = iy + incy;
		y[iy] = y[iy] + temp1 * A[i + j * lda];
		temp2 = temp2 + A[i + j * lda] * x[ix];
	    }
	    y[jy] = y[jy] + alpha * temp2;
	    jx = jx + incx;
	    jy = jy + incy;
	}
    }
    return;
}
