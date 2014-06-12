/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rorgtr.cpp,v 1.6 2009/09/22 22:46:17 nakatamaho Exp $ 
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
Rorgtr(const char *uplo, mpackint n, mpf_class * A, mpackint lda, mpf_class * tau,
    mpf_class * work, mpackint lwork, mpackint *info)
{

    mpf_class Zero = 0.0, One = 1.0;
    mpackint lquery, lwkopt, iinfo, upper, nb;
    mpackint i, j;

    *info = 0;
    if (lwork == -1)
	lquery = 1;
    else
	lquery = 0;

    upper = Mlsame_gmp(uplo, "U");
    if (!upper && !Mlsame_gmp(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((mpackint)1, n)) {
	*info = -4;
    } else {
	if (lwork < max((mpackint)1, n - 1) && !lquery) {
	    *info = -7;
	}
    }
    if (*info == 0) {
	if (upper) {
	    nb = iMlaenv_gmp(1, "Rorgql", " ", n - 1, n - 1, n - 1, -1);
	} else {
	    nb = iMlaenv_gmp(1, "Rorgqr", " ", n - 1, n - 1, n - 1, -1);
	}
	lwkopt = max((mpackint)1, n - 1) * nb;
	work[0] = (double)lwkopt;	//needs cast from double to mpf
    }
    if (*info != 0) {
	Mxerbla_gmp("Rorgtr", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	work[0] = One;
	return;
    }
    if (upper) {
//Q was determined by a call to DSYTRD with UPLO = 'U'
//Shift the vectors which define the elementary reflectors one
//column to the left, and set the last row and column of Q to
//those of the unit matrix
	for (j = 1; j <= n - 1; j++) {
	    for (i = 1; i <= j - 1; i++) {
		A[(i - 1) + (j - 1) * lda] = A[(i - 1) + j * lda];
	    }
	    A[(n - 1) + (j - 1) * lda] = Zero;
	}
	for (i = 1; i <= n - 1; i++) {
	    A[(i - 1) + (n - 1) * lda] = Zero;
	}
	A[(n - 1) + (n - 1) * lda] = One;
//Generate Q(1:n-1,1:n-1)
	Rorgql(n - 1, n - 1, n - 1, A, lda, tau, work, lwork, &iinfo);
    } else {
//Q was determined by a call to DSYTRD with UPLO = 'L'.
//Shift the vectors which define the elementary reflectors one
//column to the right, and set the first row and column of Q to
//those of the unit matrix
	for (j = n; j >= 2; j--) {
	    A[0 + (j - 1) * lda] = Zero;
	    for (i = j + 1; i <= n; i++) {
		A[(i - 1) + (j - 1) * lda] = A[(i - 1) + (j - 2) * lda];
	    }
	}
	A[0 + 0 * lda] = One;
	for (i = 2; i <= n; i++) {
	    A[(i - 1) + 0 * lda] = Zero;
	}
	if (n > 1) {
//Generate Q(2:n,2:n)
	    Rorgqr(n - 1, n - 1, n - 1, &A[1 + (1 * lda)], lda, tau,
		work, lwork, &iinfo);
	}
    }
    work[0] = (double)lwkopt;	//needs cast from double to mpf
    return;
}
