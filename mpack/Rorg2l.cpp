/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rorg2l.cpp,v 1.6 2009/09/22 22:46:17 nakatamaho Exp $ 
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
Rorg2l(mpackint m, mpackint n, mpackint k, mpf_class * A, mpackint lda, mpf_class * tau,
    mpf_class * work, mpackint *info)
{
    mpackint i, ii, j, l;
    mpf_class Zero = 0.0, One = 1.0;

    *info = 0;
    if (m < 0) {
	*info = -1;
    } else if (n < 0 || n > m) {
	*info = -2;
    } else if (k < 0 || k > n) {
	*info = -3;
    } else if (lda < max((mpackint)1, m)) {
	*info = -5;
    }
    if (*info != 0) {
	Mxerbla_gmp("Rorg2l", -(*info));
	return;
    }
//quick return if possible
    if (n <= 0)
	return;

//Initialise columns 1:n-k to columns of the unit matrix
    for (j = 0; j < n - k; j++) {
	for (l = 0; l < m; l++) {
	    A[l + j * lda] = Zero;
	}
	A[m - n + j + j * lda] = One;
    }

    for (i = 1; i <= k; i++) {
	ii = n - k + i;
//Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
	A[(m - n + ii - 1) + (ii - 1) * lda] = One;
	Rlarf("Left", m - n + ii, ii - 1, &A[0 + (ii - 1) * lda], 1,
	    tau[i - 1], A, lda, work);
	Rscal(m - n + ii - 1, -tau[i - 1], &A[0 + (ii - 1) * lda], 1);
	A[(m - n + ii - 1) + (ii - 1) * lda] = One - tau[i - 1];
//Set A(m-k+i+1:m,n-k+i) to zero
	for (l = m - n + ii + 1; l <= m; l++) {
	    A[(l - 1) + (ii - 1) * lda] = Zero;
	}
    }
    return;
}
