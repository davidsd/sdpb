/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rsyr2.cpp,v 1.4 2009/09/24 07:25:57 nakatamaho Exp $ 
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
Based on http://www.netlib.org/blas/dsyr2.f
Rsyr2 performs the symmetric rank 2 operation
 A := alpha*x*y' + alpha*y*x' + A,
where alpha is a scalar, x and y are n element vectors and A is an n
by n symmetric matrix.
*/

#include <mblas_gmp.h>

void
Rsyr2(const char *uplo, mpackint n, mpf_class alpha, mpf_class * x,
    mpackint incx, mpf_class * y, mpackint incy, mpf_class * A, mpackint lda)
{
    mpf_class temp1, temp2;

    mpf_class Zero = 0.0;

    mpackint ix, iy, jx, jy, kx, ky;

    //test the input parameters.
    mpackint info = 0;

    if (!Mlsame_gmp(uplo, "U") && !Mlsame_gmp(uplo, "L"))
	info = 1;
    else if (n < 0)
	info = 2;
    else if (incx == 0)
	info = 5;
    else if (incy == 0)
	info = 7;
    else if (lda < max((mpackint) 1, n))
	info = 9;
    if (info != 0) {
	Mxerbla_gmp("Rsyr2 ", info);
	return;
    }
    //quick return if possible.
    if ((n == 0) || (alpha == Zero))
	return;

    if (incx > 0)
	kx = 0;
    else
	kx = -(n - 1) * incx;
    if (incy > 0)
	ky = 0;
    else
	ky = -(n - 1) * incy;
    jx = kx;
    jy = ky;

    if (Mlsame_gmp(uplo, "U")) {
	for (mpackint j = 0; j < n; j++) {
	    if ((x[jx] != Zero) || (y[jy] != Zero)) {
		temp1 = alpha * y[jy];
		temp2 = alpha * x[jx];
		ix = kx;
		iy = ky;
		for (mpackint i = 0; i <= j; i++) {
		    A[i + j * lda] =
			A[i + j * lda] + x[ix] * temp1 + y[iy] * temp2;
		    ix = ix + incx;
		    iy = iy + incy;
		}
	    }
	    jx = jx + incx;
	    jy = jy + incy;
	}
    } else {
	//form  a  when a is stored in the lower triangle.
	for (mpackint j = 0; j < n; j++) {
	    if ((x[jx] != Zero) || (y[jy] != Zero)) {
		temp1 = alpha * y[jy];
		temp2 = alpha * x[jx];
		ix = jx;
		iy = jy;
		for (mpackint i = j; i < n; i++) {
		    A[i + j * lda] =
			A[i + j * lda] + x[ix] * temp1 + y[iy] * temp2;
		    ix = ix + incx;
		    iy = iy + incy;
		}
	    }
	    jx = jx + incx;
	    jy = jy + incy;
	}
    }
    return;
}
