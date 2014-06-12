/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rlansy.cpp,v 1.2 2009/09/12 07:59:10 nakatamaho Exp $ 
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

mpf_class
Rlansy(const char *norm, const char *uplo, mpackint n, mpf_class * A, mpackint lda,
    mpf_class * work)
{
    mpf_class One = 1.0, Zero = 0.0;
    mpackint i, j;
    mpf_class absa, scale, sum, value;
    mpf_class mtmp;

    if (n == 0) {
	value = Zero;
	return value;
    }
    if (Mlsame_gmp(norm, "M")) {
//Find max(abs(A(i,j))).
	value = Zero;
	if (Mlsame_gmp(uplo, "U")) {
	    for (j = 0; j < n; j++) {
		for (i = 0; i <= j; i++) {
		    mtmp = abs(A[i + j * lda]);
		    value = max(value, mtmp);
		}
	    }
	} else {
	    for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
		    mtmp = abs(A[i + j * lda]);
		    value = max(value, mtmp);
		}
	    }
	}
    } else if (Mlsame_gmp(norm, "I") || Mlsame_gmp(norm, "O") || Mlsame_gmp(norm, "1")) {
// Find normI(A) ( = norm1(A), since A is symmetric).
	value = Zero;
	if (Mlsame_gmp(uplo, "U")) {
	    for (j = 0; j < n; j++) {
		sum = Zero;
		for (i = 0; i < j; i++) {
		    absa = abs(A[i + j * lda]);
		    sum += absa;
		    work[i] += absa;
		}
		work[j] = sum + abs(A[j + j * lda]);
	    }
	    for (i = 0; i < n; i++) {
		value = max(value, work[i]);
	    }
	} else {
	    for (i = 0; i < n; i++) {
		work[i] = Zero;
	    }
	    for (j = 0; j < n; j++) {
		sum = work[j] + abs(A[j + j * lda]);
		for (i = j + 1; i < n; i++) {
		    absa = abs(A[i + j * lda]);
		    sum += absa;
		    work[i] += absa;
		}
		value = max(value, sum);
	    }
	}
    } else if (Mlsame_gmp(norm, "F") || Mlsame_gmp(norm, "E")) {
//Find normF(A).
	scale = Zero;
	sum = One;
	if (Mlsame_gmp(uplo, "U")) {
	    for (j = 1; j < n; j++) {
		Rlassq(j, &A[j * lda], 1, &scale, &sum);
	    }
	} else {
	    for (j = 0; j < n - 1; j++) {
		Rlassq(n - j - 1, &A[(j + 1) + j * lda], 1, &scale, &sum);
	    }
	}
	sum *= 2.0;
	Rlassq(n, A, lda + 1, &scale, &sum);
	value = scale * sqrt(sum);
    }
    return value;
}
