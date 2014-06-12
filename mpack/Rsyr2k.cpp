/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rsyr2k.cpp,v 1.4 2009/09/24 07:25:57 nakatamaho Exp $ 
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
http://www.netlib.org/blas/dsyr2k.f
Rsyr2k performs one of the symmetric rank 2k operations
C := alpha*A*B' + alpha*B*A' + beta*C,
 or
C := alpha*A'*B + alpha*B'*A + beta*C,
where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
and  A and B  are  n by k  matrices  in the  first  case  and  k by n
matrices in the second case.
*/

#include <mblas_gmp.h>

void
Rsyr2k(const char *uplo, const char *trans, mpackint n, mpackint k,
    mpf_class alpha, mpf_class * A, mpackint lda, mpf_class * B, mpackint ldb,
    mpf_class beta, mpf_class * C, mpackint ldc)
{
    mpackint nrowa, upper, info;

    mpf_class Zero = 0.0, One = 1.0;

    mpf_class temp1, temp2;

    //test the input parameters.
    if (Mlsame_gmp(trans, "N"))
	nrowa = n;
    else
	nrowa = k;
    upper = Mlsame_gmp(uplo, "U");

    info = 0;
    if ((!upper) && (!Mlsame_gmp(uplo, "L")))
	info = 1;
    else if ((!Mlsame_gmp(trans, "N")) && (!Mlsame_gmp(trans, "T"))
	&& (!Mlsame_gmp(trans, "C")))
	info = 2;
    else if (n < 0)
	info = 3;
    else if (k < 0)
	info = 4;
    else if (lda < max((mpackint) 1, nrowa))
	info = 7;
    else if (ldb < max((mpackint) 1, nrowa))
	info = 9;
    else if (ldc < max((mpackint) 1, n))
	info = 12;
    if (info != 0) {
	Mxerbla_gmp("Rsyr2k", info);
	return;
    }
    //quick return if possible.
    if ((n == 0) || (((alpha == Zero) || (k == 0)) && (beta == One)))
	return;

    //and when alpha==Zero.
    if (alpha == Zero) {
	if (upper) {
	    if (beta == Zero) {
		for (mpackint j = 0; j < n; j++) {
		    for (mpackint i = 0; i <= j; i++) {
			C[i + j * ldc] = Zero;
		    }
		}
	    } else {
		for (mpackint j = 0; j < n; j++) {
		    for (mpackint i = 0; i <= j; i++) {
			C[i + j * ldc] = beta * C[i + j * ldc];
		    }
		}
	    }
	} else {
	    if (beta == Zero) {
		for (mpackint j = 0; j < n; j++) {
		    for (mpackint i = j; i < n; i++) {
			C[i + j * ldc] = Zero;
		    }
		}
	    } else {
		for (mpackint j = 0; j < n; j++) {
		    for (mpackint i = j; i < n; i++) {
			C[i + j * ldc] = beta * C[i + j * ldc];
		    }
		}
	    }
	}
	return;
    }
    //start the operations.
    if (Mlsame_gmp(trans, "N")) {
	//form C:= alpha*A*B' + alpha*B*A'+C.
	if (upper) {
	    for (mpackint j = 0; j < n; j++) {
		if (beta == Zero) {
		    for (mpackint i = 0; i <= j; i++) {
			C[i + j * ldc] = Zero;
		    }
		} else if (beta != One) {
		    for (mpackint i = 0; i <= j; i++) {
			C[i + j * ldc] = beta * C[i + j * ldc];
		    }
		}
		for (mpackint l = 0; l < k; l++) {
		    if ((A[j + l * lda] != Zero) || (B[j + l * ldb] != Zero)) {
			temp1 = alpha * B[j + l * ldb];
			temp2 = alpha * A[j + l * lda];
			for (mpackint i = 0; i <= j; i++) {
			    C[i + j * ldc] =
				C[i + j * ldc] + A[i + l * lda] * temp1 + B[i +
				l * ldb] * temp2;
			}
		    }
		}
	    }
	} else {
	    for (mpackint j = 0; j < n; j++) {
		if (beta == Zero) {
		    for (mpackint i = j; i < n; i++) {
			C[i + j * ldc] = Zero;
		    }
		} else if (beta != One) {
		    for (mpackint i = j; i < n; i++) {
			C[i + j * ldc] = beta * C[i + j * ldc];
		    }
		}
		for (mpackint l = 0; l < k; l++) {
		    if ((A[j + l * lda] != Zero) || (B[j + l * ldb] != Zero)) {
			temp1 = alpha * B[j + l * ldb];
			temp2 = alpha * A[j + l * lda];
			for (mpackint i = j; i < n; i++) {
			    C[i + j * ldc] =
				C[i + j * ldc] + A[i + l * lda] * temp1 + B[i +
				l * ldb] * temp2;
			}
		    }
		}
	    }
	}
    } else {
	//form  C := alpha*A'*B + alpha*B'*A + C.
	if (upper) {
	    for (mpackint j = 0; j < n; j++) {
		for (mpackint i = 0; i <= j; i++) {
		    temp1 = Zero;
		    temp2 = Zero;
		    for (mpackint l = 0; l < k; l++) {
			temp1 = temp1 + A[l + i * lda] * B[l + j * ldb];
			temp2 = temp2 + B[l + i * ldb] * A[l + j * lda];
		    }
		    if (beta == Zero) {
			C[i + j * ldc] = alpha * temp1 + alpha * temp2;
		    } else {
			C[i + j * ldc] =
			    beta * C[i + j * ldc] + alpha * temp1 +
			    alpha * temp2;
		    }
		}
	    }
	} else {
	    for (mpackint j = 0; j < n; j++) {
		for (mpackint i = j; i < n; i++) {
		    temp1 = Zero;
		    temp2 = Zero;
		    for (mpackint l = 0; l < k; l++) {
			temp1 = temp1 + A[l + i * lda] * B[l + j * ldb];
			temp2 = temp2 + B[l + i * ldb] * A[l + j * lda];
		    }
		    if (beta == Zero) {
			C[i + j * ldc] = alpha * temp1 + alpha * temp2;
		    } else {
			C[i + j * ldc] =
			    beta * C[i + j * ldc] + alpha * temp1 +
			    alpha * temp2;
		    }
		}
	    }
	}
    }
    return;
}
