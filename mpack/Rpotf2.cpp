/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rpotf2.cpp,v 1.7 2009/09/25 04:00:39 nakatamaho Exp $ 
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
#include <stdlib.h>

void
Rpotf2(const char *uplo, mpackint n, mpf_class * A, mpackint lda, mpackint *info)
{
    mpackint j, upper, success = 1;
    mpf_class ajj;
    mpf_class Zero = 0.0;
    mpf_class One = 1.0;

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
	Mxerbla_gmp("Rpotf2", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0)
	return;

    if (upper) {
//Compute the Cholesky factorization A = U'*U.
	for (j = 0; j < n; j++) {
//Compute U(J,J) and test for non-positive-definiteness.
	    ajj = A[j + j * lda] - Rdot(j, &A[j * lda], 1, &A[j * lda], 1);
	    if (ajj <= Zero) {
		A[j + j * lda] = ajj;
		success = 0;
		break;
	    }
	    ajj = sqrt(ajj);
	    A[j + j * lda] = ajj;
//Compute elements J+1:N of row J.
	    if (j < n) {
		Rgemv("Transpose", j, n - j - 1, -One, &A[(j + 1) * lda], lda,
		    &A[j * lda], 1, One, &A[j + (j + 1) * lda], lda);
		Rscal(n - j - 1, One / ajj, &A[j + (j + 1) * lda], lda);
	    }
	}
    } else {
//Compute the Cholesky factorization A = L*L'.
	for (j = 0; j < n; j++) {
// Compute L(J,J) and test for non-positive-definiteness.
	    ajj = A[j + j * lda] - Rdot(j, &A[j], lda, &A[j], lda);
	    if (ajj <= Zero) {
		A[j + j * lda] = ajj;
		success = 0;
		break;
	    }
	    ajj = sqrt(ajj);
	    A[j + j * lda] = ajj;

//Compute elements J+1:N of column J.
	    if (j < n) {
		Rgemv("No transpose", n - j - 1, j, -One, &A[j + 1], lda,
		    &A[j], lda, One, &A[j + 1 + j * lda], 1);
		Rscal(n - j - 1, One / ajj, &A[j + 1 + j * lda], 1);
	    }
	}
    }
    if (!success)
	*info = j + 1;
    return;
}
