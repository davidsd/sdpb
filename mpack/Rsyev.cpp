/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rsyev.cpp,v 1.7 2009/09/22 21:22:09 nakatamaho Exp $ 
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
Rsyev(const char *jobz, const char *uplo, mpackint n, mpf_class * A,
    mpackint lda, mpf_class * w, mpf_class * work, mpackint *lwork, mpackint *info)
{

    mpackint wantz, lower, lquery, nb, lwkopt, iscale, imax;
    mpackint inde, indtau, indwrk, llwork, iinfo;

    mpf_class Zero = 0.0, One = 1.0, Two = 2.0;
    mpf_class safmin, eps, smlnum, bignum, rmin, rmax;
    mpf_class sigma, anrm;
    mpf_class rtmp;

    wantz = Mlsame_gmp(jobz, "V");
    lower = Mlsame_gmp(uplo, "L");
    lquery = 0;
    if (*lwork == -1)
	lquery = 1;

    *info = 0;
    if (!(wantz || Mlsame_gmp(jobz, "N"))) {
	*info = -1;
    } else if (!(lower || Mlsame_gmp(uplo, "U"))) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (lda < max((mpackint)1, n)) {
	*info = -5;
    }

    if (*info == 0) {
	nb = iMlaenv_gmp(1, "Rsytrd", uplo, n, -1, -1, -1);
	lwkopt = max((mpackint)1, (nb + 2) * n);
	work[0] = (double)lwkopt;	//needs cast mpackint to mpf
	if (*lwork < max((mpackint)1, 3 * n - 1) && !lquery) {
	    *info = -8;
	}
    }

    if (*info != 0) {
	Mxerbla_gmp("Rsyev ", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
    if (n == 1) {
	w[0] = A[0];
	work[0] = Two;
	if (wantz) {
	    A[0] = One;
	}
	return;
    }
//Get machine constants.
    safmin = Rlamch_gmp("Safe minimum");
    eps = Rlamch_gmp("Precision");
    smlnum = safmin / eps;
    bignum = One / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);

//Scale matrix to allowable range, if necessary.
    anrm = Rlansy("M", uplo, n, A, lda, work);
    iscale = 0;
    if (anrm > Zero && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	Rlascl(uplo, 0, 0, One, sigma, n, n, A, lda, info);
    }
//Call DSYTRD to reduce symmetric matrix to tridiagonal form.
    inde = 1;
    indtau = inde + n;
    indwrk = indtau + n;
    llwork = *lwork - indwrk + 1;
    Rsytrd(uplo, n, &A[0], lda, &w[0], &work[inde - 1], &work[indtau - 1],
	&work[indwrk - 1], llwork, &iinfo);

//For eigenvalues only, call DSTERF.  For eigenvectors, first call
//DORGTR to generate the orthogonal matrix, then call DSTEQR.
    if (!wantz) {
	Rsterf(n, &w[0], &work[inde - 1], info);
    } else {
	Rorgtr(uplo, n, A, lda, &work[indtau - 1], &work[indwrk - 1], llwork,
	    &iinfo);
	Rsteqr(jobz, n, w, &work[inde - 1], A, lda, &work[indtau - 1], info);
    }

//If matrix was scaled, then rescale eigenvalues appropriately.
    if (iscale == 1) {
	if (*info == 0) {
	    imax = n;
	} else {
	    imax = *info - 1;
	}
	rtmp = One / sigma;
	Rscal(imax, rtmp, &w[0], 1);
    }
//Set WORK(1) to optimal workspace size.
    work[0] = (double)lwkopt;	//needs cast from mpackint to mpf

    return;
}
