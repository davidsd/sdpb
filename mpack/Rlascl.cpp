/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rlascl.cpp,v 1.7 2009/09/25 04:00:39 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void
Rlascl(const char *type, mpackint kl, mpackint ku, mpf_class cfrom, mpf_class cto, mpackint m,
    mpackint n, mpf_class * A, mpackint lda, mpackint *info)
{
    mpackint i, j, k1, k2, k3, k4;
    mpackint itype;
    mpf_class One = 1.0, Zero = 0.0;
    mpf_class bignum, cfrom1, cfromc, cto1, ctoc, mul, smlnum;
    mpackint done = MFALSE;

    *info = 0;
    if (Mlsame_gmp(type, "G")) {
	itype = 0;
    } else if (Mlsame_gmp(type, "L")) {
	itype = 1;
    } else if (Mlsame_gmp(type, "U")) {
	itype = 2;
    } else if (Mlsame_gmp(type, "H")) {
	itype = 3;
    } else if (Mlsame_gmp(type, "B")) {
	itype = 4;
    } else if (Mlsame_gmp(type, "Q")) {
	itype = 5;
    } else if (Mlsame_gmp(type, "Z")) {
	itype = 6;
    } else {
	itype = -1;
    }
    if (itype == -1) {
	*info = -1;
    } else if (cfrom == Zero) {
	*info = -4;
    } else if (m < 0) {
	*info = -6;
    } else if (n < 0 || (itype == 4 && n != m) || (itype == 5 && n != m) ) {
	*info = -7;
    } else if (itype <= 3 && lda < max((mpackint)1, m)) {
	*info = -9;
    } else if (itype >= 4) {
	if (kl < 0 || kl > max(m - 1, (mpackint)0)) {
	    *info = -2;
	} else {
	    if (ku < 0 || ku > max(n - 1, (mpackint)0) || ((itype == 4 || itype == 5) &&
		kl != ku)) {
		*info = -3;
	    } else if ( (itype == 4 && lda < kl + 1) || (itype == 5 && lda < ku + 1)
		|| (itype == 6 && lda < (kl * 2) + ku + 1)) {
		*info = -9;
	    }
	}
    }

    if (*info != 0) {
	Mxerbla_gmp("Rlascl", -(*info));
	return;
    }
//Quick return if possible 
    if (n == 0 || m == 0) {
	return;
    }
//Get machine parameters 
    smlnum = Rlamch_gmp("S");
    bignum = One / smlnum;

    cfromc = cfrom;
    ctoc = cto;

    while (done == MFALSE) {
	cfrom1 = cfromc * smlnum;
	cto1 = ctoc / bignum;
	if (abs(cfrom1) > abs(ctoc) && ctoc != Zero) {
	    mul = smlnum;
	    done = MFALSE;
	    cfromc = cfrom1;
	} else if (abs(cto1) > abs(cfromc)) {
	    mul = bignum;
	    done = MFALSE;
	    ctoc = cto1;
	} else {
	    mul = ctoc / cfromc;
	    done = MTRUE;
	}
	if (itype == 0) {
//Full matrix
	    for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
		    A[i + j * lda] = A[i + j * lda] * mul;
		}
	    }
	} else if (itype == 1) {
//Lower triangular matrix
	    for (j = 0; j < n; j++) {
		for (i = j; i < m; i++) {
		    A[i + j * lda] = A[i + j * lda] * mul;
		}
	    }
	} else if (itype == 2) {
//Upper triangular matrix
	    for (j = 0; j < n; j++) {
		for (i = 0; i <= min(j, m - 1); i++) {
		    A[i + j * lda] = A[i + j * lda] * mul;
		}
	    }
	} else if (itype == 3) {
//Upper Hessenberg matrix
	    for (j = 0; j < n; j++) {
		for (i = 0; i <= min(j + 1, m - 1); i++) {
		    A[i + j * lda] = A[i + j * lda] * mul;
		}
	    }
	} else if (itype == 4) {
//Lower half of a symmetric band matrix
	    k3 = kl + 1;
	    k4 = n + 1;
	    for (j = 0; j < n; j++) {
		for (i = 0; i < min(k3, k4 - j - 1); i++) {
		    A[i + j * lda] *= mul;
		}
	    }

	} else if (itype == 5) {
//Upper half of a symmetric band matrix
	    k1 = ku + 2;
	    k3 = ku + 1;
	    for (j = 0; j < n; j++) {
		for (i = max(k1 - j - 1, (mpackint)1) - 1; i < k3; i++) {
		    A[i + j * lda] = A[i + j * lda] * mul;
		}
	    }
	} else if (itype == 6) {
//Band matrix
	    k1 = kl + ku + 2;
	    k2 = kl + 1;
	    k3 = (kl << 1) + ku + 1;
	    k4 = kl + ku + 1 + m;
	    for (j = 0; j < n; j++) {
		for (i = max(k1 - j - 1, k2) - 1; i < min(k3, k4 - j - 1); i++) {
		    A[i + j * lda] = A[i + j * lda] * mul;
		}
	    }
	}
    }
    return;
}
