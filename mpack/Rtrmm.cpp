/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rtrmm.cpp,v 1.4 2009/09/24 07:25:57 nakatamaho Exp $ 
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
Based on http://www.netlib.org/blas/dtrmm.f
Rtrmm performs one of the matrix-matrix operations
 B := alpha*op(A)*B, or B := alpha*B*op(A),
where alpha is a scalar, B is an m by n matrix, A is a unit, or
non-unit, upper or lower triangular matrix and op(A) is one  of
 op(A) = A  or op(A) = A'.
*/

#include <mblas_gmp.h>

void
Rtrmm(const char *side, const char *uplo, const char *transa, const char *diag,
    mpackint m, mpackint n, mpf_class alpha, mpf_class * A, mpackint lda,
    mpf_class * B, mpackint ldb)
{
    mpackint info, lside, nrowa, nounit, upper;

    mpf_class temp;

    mpf_class Zero = 0.0, One = 1.0;

    //test the input parameters.
    lside = Mlsame_gmp(side, "L");
    if (lside)
	nrowa = m;
    else
	nrowa = n;

    nounit = Mlsame_gmp(diag, "N");
    upper = Mlsame_gmp(uplo, "U");

    info = 0;
    if ((!lside) && (!Mlsame_gmp(side, "R")))
	info = 1;
    else if ((!upper) && (!Mlsame_gmp(uplo, "L")))
	info = 2;
    else if ((!Mlsame_gmp(transa, "N")) && (!Mlsame_gmp(transa, "T"))
	&& (!Mlsame_gmp(transa, "C")))
	info = 3;
    else if ((!Mlsame_gmp(diag, "U")) && (!Mlsame_gmp(diag, "N")))
	info = 4;
    else if (m < 0)
	info = 5;
    else if (n < 0)
	info = 6;
    else if (lda < max((mpackint) 1, nrowa))
	info = 9;
    else if (ldb < max((mpackint) 1, m))
	info = 11;
    if (info != 0) {
	Mxerbla_gmp("Rtrmm ", info);
	return;
    }
    //quick return if possible.
    if (m == 0 || n == 0)
	return;

    //and when alpha==Zero.
    if (alpha == Zero) {
	for (mpackint j = 0; j < n; j++) {
	    for (mpackint i = 0; i < m; i++) {
		B[i + j * ldb] = Zero;
	    }
	}
	return;
    }
    //start the operations.
    if (lside) {
	if (Mlsame_gmp(transa, "N")) {
	    //Form B := alpha*A*B.
	    if (upper) {
		for (mpackint j = 0; j < n; j++) {
		    for (mpackint k = 0; k < m; k++) {
			if (B[k + j * ldb] != Zero) {
			    temp = alpha * B[k + j * ldb];
			    for (mpackint i = 0; i < k; i++) {
				B[i + j * ldb] =
				    B[i + j * ldb] + temp * A[i + k * lda];
			    }
			    if (nounit)
				temp = temp * A[k + k * lda];
			    B[k + j * ldb] = temp;
			}
		    }
		}
	    } else {
		for (mpackint j = 0; j < n; j++) {
		    for (mpackint k = m - 1; k >= 0; k--) {
			if (B[k + j * ldb] != Zero) {
			    temp = alpha * B[k + j * ldb];
			    B[k + j * ldb] = temp;
			    if (nounit)
				B[k + j * ldb] =
				    B[k + j * ldb] * A[k + k * lda];
			    for (mpackint i = k + 1; i < m; i++) {
				B[i + j * ldb] =
				    B[i + j * ldb] + temp * A[i + k * lda];
			    }
			}
		    }
		}
	    }
	} else {
	    //Form B := alpha*A'*B.
	    if (upper) {
		for (mpackint j = 0; j < n; j++) {
		    for (mpackint i = m - 1; i >= 0; i--) {
			temp = B[i + j * ldb];
			if (nounit)
			    temp = temp * A[i + i * lda];
			for (mpackint k = 0; k < i; k++) {
			    temp = temp + A[k + i * lda] * B[k + j * ldb];
			}
			B[i + j * ldb] = alpha * temp;
		    }
		}
	    } else {
		for (mpackint j = 0; j < n; j++) {
		    for (mpackint i = 0; i < m; i++) {
			temp = B[i + j * ldb];
			if (nounit)
			    temp = temp * A[i + i * lda];
			for (mpackint k = i + 1; k < m; k++) {
			    temp = temp + A[k + i * lda] * B[k + j * ldb];
			}
			B[i + j * ldb] = alpha * temp;
		    }
		}
	    }
	}
    } else {
	if (Mlsame_gmp(transa, "N")) {
	    //Form B := alpha*B*A.
	    if (upper) {
		for (mpackint j = n - 1; j >= 0; j--) {
		    temp = alpha;
		    if (nounit)
			temp = temp * A[j + j * lda];
		    for (mpackint i = 0; i < m; i++) {
			B[i + j * ldb] = temp * B[i + j * ldb];
		    }
		    for (mpackint k = 0; k < j; k++) {
			if (A[k + j * lda] != Zero) {
			    temp = alpha * A[k + j * lda];
			    for (mpackint i = 0; i < m; i++) {
				B[i + j * ldb] =
				    B[i + j * ldb] + temp * B[i + k * ldb];
			    }
			}
		    }
		}
	    } else {
		for (mpackint j = 0; j < n; j++) {
		    temp = alpha;
		    if (nounit)
			temp = temp * A[j + j * lda];
		    for (mpackint i = 0; i < m; i++) {
			B[i + j * ldb] = temp * B[i + j * ldb];
		    }
		    for (mpackint k = j + 1; k < n; k++) {
			if (A[k + j * lda] != Zero) {
			    temp = alpha * A[k + j * lda];
			    for (mpackint i = 0; i < m; i++) {
				B[i + j * ldb] =
				    B[i + j * ldb] + temp * B[i + k * ldb];
			    }
			}
		    }
		}
	    }
	} else {
	    if (upper) {
		for (mpackint k = 0; k < n; k++) {
		    for (mpackint j = 0; j < k; j++) {
			if (A[j + k * lda] != Zero) {
			    temp = alpha * A[j + k * lda];
			    for (mpackint i = 0; i < m; i++) {
				B[i + j * ldb] =
				    B[i + j * ldb] + temp * B[i + k * ldb];
			    }
			}
		    }
		    temp = alpha;
		    if (nounit)
			temp = temp * A[k + k * lda];
		    if (temp != One) {
			for (mpackint i = 0; i < m; i++) {
			    B[i + k * ldb] = temp * B[i + k * ldb];
			}
		    }
		}
	    } else {
		for (mpackint k = n - 1; k >= 0; k--) {
		    for (mpackint j = k + 1; j < n; j++) {
			if (A[j + k * lda] != Zero) {
			    temp = alpha * A[j + k * lda];
			    for (mpackint i = 0; i < m; i++) {
				B[i + j * ldb] =
				    B[i + j * ldb] + temp * B[i + k * ldb];
			    }
			}
		    }
		    temp = alpha;
		    if (nounit)
			temp = temp * A[k + k * lda];
		    if (temp != One) {
			for (mpackint i = 0; i < m; i++) {
			    B[i + k * ldb] = temp * B[i + k * ldb];
			}
		    }
		}
	    }
	}
    }
    return;
}
