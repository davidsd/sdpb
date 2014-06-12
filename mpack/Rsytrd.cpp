/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rsytrd.cpp,v 1.7 2009/09/22 21:22:09 nakatamaho Exp $ 
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
Rsytrd(const char *uplo, mpackint n, mpf_class * A, mpackint lda, mpf_class * d,
    mpf_class * e, mpf_class * tau, mpf_class * work, mpackint lwork, mpackint *info)
{
    mpackint upper, lquery, nb, lwkopt, nx, iws;
    mpackint ldwork, nbmin, kk;
    mpackint i, j;
    mpackint iinfo;
    mpf_class One = 1.0;

    *info = 0;
    upper = Mlsame_gmp(uplo, "U");
    lquery = 0;
    if (lwork == -1)
	lquery = 1;

    if (!upper && !Mlsame_gmp(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((mpackint)1, n)) {
	*info = -4;
    } else if (lwork < 1 && !lquery) {
	*info = -9;
    }
    if (*info == 0) {
//Determine the block size.
	nb = iMlaenv_gmp(1, "Rsytrd", uplo, n, -1, -1, -1);
	lwkopt = n * nb;
	work[0] = (double)lwkopt;	//cast from mpackint to mpf
    }
    if (*info != 0) {
	Mxerbla_gmp("Rsytrd", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	work[0] = One;
	return;
    }

    nx = n;
    iws = 1;
    if (nb > 1 && nb < n) {
//Determine when to cross over from blocked to unblocked code
//(last block is always handled by unblocked code).
	nx = max(nb, iMlaenv_gmp(3, "Rsytrd", uplo, n, -1, -1, -1));
	if (nx < n) {
//Determine if workspace is large enough for blocked code.
	    ldwork = n;
	    iws = ldwork * nb;
	    if (lwork < iws) {
//Not enough workspace to use optimal NB:  determine the
//minimum value of NB, and reduce NB or force use of
//unblocked code by setting NX = N.
		nb = max(lwork / ldwork, (mpackint)1);
		nbmin = iMlaenv_gmp(2, "Rsytrd", uplo, n, -1, -1, -1);
		if (nb < nbmin) {
		    nx = n;
		}
	    }
	} else {
	    nx = n;
	}
    } else {
	nb = 1;
    }
    if (upper) {
//Reduce the upper triangle of A.
//Columns 1:kk are handled by the unblocked method.
	kk = n - ((n - nx + nb - 1) / nb) * nb;
	for (i = n - nb + 1; i >= kk + 1; i = i - nb) {
// Reduce columns i:i+nb-1 to tridiagonal form and form the
//matrix W which is needed to update the unreduced part of
//the matrix
	    Rlatrd(uplo, i + nb - 1, nb, A, lda, e, tau, work, ldwork);
//Update the unreduced submatrix A(1:i-1,1:i-1), using an
//update of the form:  A := A - V*W' - W*V'
	    Rsyr2k(uplo, "No transpose", i - 1, nb, -One,
		&A[0 + (i - 1) * lda], lda, work, ldwork, One, A, lda);
//Copy superdiagonal elements back into A, and diagonal
//elements into D
	    for (j = i; j <= i + nb - 1; j++) {
		A[(j - 2) + (j - 1) * lda] = e[j - 2];
		d[j - 1] = A[(j - 1) + (j - 1) * lda];
	    }
	}
//Use unblocked code to reduce the last or only block
	Rsytd2(uplo, kk, A, lda, d, e, tau, &iinfo);
    } else {
//Reduce the lower triangle of A
	for (i = 1; i <= n - nx; i = i + nb) {
//Reduce columns i:i+nb-1 to tridiagonal form and form the
//matrix W which is needed to update the unreduced part of
//the matrix
	    Rlatrd(uplo, n - i + 1, nb, &A[(i - 1) + (i - 1) * lda], lda,
		&e[i - 1], &tau[i - 1], work, ldwork);
//Update the unreduced submatrix A(i+ib:n,i+ib:n), using
//an update of the form:  A := A - V*W' - W*V'
	    Rsyr2k(uplo, "No transpose", n - i - nb + 1, nb, -One,
		&A[(i + nb - 1) + (i - 1) * lda], lda, &work[nb], ldwork, One,
		&A[(i + nb - 1) + (i + nb - 1) * lda], lda);
//Copy subdiagonal elements back into A, and diagonal
//elements into D
	    for (j = i; j <= i + nb - 1; j++) {
		A[j + (j - 1) * lda] = e[j - 1];
		d[j - 1] = A[(j - 1) + (j - 1) * lda];
	    }
	}
//Use unblocked code to reduce the last or only block
	Rsytd2(uplo, n - i + 1, &A[(i - 1) + (i - 1) * lda], lda, &d[i - 1],
	    &e[i - 1], &tau[i - 1], &iinfo);
    }
    work[0] = (double)lwkopt;	//cast mpf to mpackint
    return;
}
