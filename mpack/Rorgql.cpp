/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rorgql.cpp,v 1.7 2009/09/25 04:00:39 nakatamaho Exp $ 
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
Rorgql(mpackint m, mpackint n, mpackint k, mpf_class * A, mpackint lda, mpf_class * tau,
    mpf_class * work, mpackint lwork, mpackint *info)
{
    mpf_class Zero = 0.0, One = 1.0;
    mpackint nbmin, nx, iws, nb, lwkopt, lquery, kk;
    mpackint i, j, l, iinfo, ldwork, ib;

//Test the input arguments
    *info = 0;
    if (lwork == -1)
	lquery = 1;
    else
	lquery = 0;

    if (m < 0) {
	*info = -1;
    } else if (n < 0 || n > m) {
	*info = -2;
    } else if (k < 0 || k > n) {
	*info = -3;
    } else if (lda < max((mpackint)1, m)) {
	*info = -5;
    }

    if (*info == 0) {
	if (n == 0) {
	    lwkopt = 1;
	} else {
	    nb = iMlaenv_gmp(1, "Rorgql", " ", m, n, k, -1);
	    lwkopt = n * nb;
	}
	work[0] = (double)lwkopt;	//needs cast mpackint to mpf
	if (lwork < max((mpackint)1, n) && !lquery) {
	    *info = -8;
	}
    }
    if (*info != 0) {
	Mxerbla_gmp("Rorgql", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n <= 0)
	return;
    nbmin = 2;
    nx = 0;
    iws = n;
    if (nb > 1 && nb < k) {
//Determine when to cross over from blocked to unblocked code.
	nx = max((mpackint)0, iMlaenv_gmp(3, "Rorgql", " ", m, n, k, -1));
	if (nx < k) {
//Determine if workspace is large enough for blocked code.
	    ldwork = n;
	    iws = ldwork * nb;
	    if (lwork < iws) {
//Not enough workspace to use optimal NB:  reduce NB and
//determine the minimum value of NB.
		nb = lwork / ldwork;
		nbmin = max((mpackint)2, iMlaenv_gmp(2, "Rorgql", " ", m, n, k, -1));
	    }
	}
    }
    if (nb >= nbmin && nb < k && nx < k) {
//Use blocked code after the first block.
//The last kk columns are handled by the block method.
	kk = min(k, (k - nx + nb - 1) / nb * nb);
//Set A(m-kk+1:m,1:n-kk) to zero.
	for (j = 1; j <= n - kk; j++) {
	    for (i = m - kk + 1; i <= m; i++) {
		A[(i - 1) + (j - 1) * lda] = Zero;
	    }
	}
    } else {
	kk = 0;
    }
//Use unblocked code for the first or only block.
    Rorg2l(m - kk, n - kk, k - kk, A, lda, tau, work, &iinfo);
    if (kk > 0) {
	for (i = k - kk + 1; i <= k; i = i + nb) {
	    ib = min(nb, k - i + 1);
	    if (n - k + i > 1) {
//Form the triangular factor of the block reflector
//H = H(i+ib-1) . . . H(i+1) H(i)
		Rlarft("Backward", "Columnwise", m - k + i + ib - 1, ib,
		    &A[0 + (n - k + i - 1) * lda], lda, &tau[i - 1], work,
		    ldwork);
//Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
		Rlarfb("Left", "No transpose", "Backward", "Columnwise",
		    m - k + i + ib - 1, n - k + i - 1, ib,
		    &A[0 + (n - k + i - 1) * lda], lda, work, ldwork, A,
		    lda, &work[ib], ldwork);
	    }
//Apply H to rows 1:m-k+i+ib-1 of current block
	    Rorg2l(m - k + i + ib - 1, ib, ib, &A[0 + (n - k + i - 1) * lda],
		lda, &tau[i - 1], work, &iinfo);
//Set rows m-k+i+ib:m of current block to zero
	    for (j = n - k + i; j <= n - k + i + ib - 1; j++) {
		for (l = m - k + i + ib; l <= m; l++) {
		    A[(l - 1) + (j - 1) * lda] = Zero;
		}
	    }
	}
    }
    work[0] = (double)iws;	//needs cast mpackint to mpf
    return;
}
