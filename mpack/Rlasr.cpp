/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rlasr.cpp,v 1.5 2009/09/22 20:33:23 nakatamaho Exp $ 
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
Rlasr(const char *side, const char *pivot, const char *direct, mpackint m,
    mpackint n, mpf_class * c, mpf_class * s, mpf_class * A, mpackint lda)
{
    mpf_class Zero = 0.0;
    mpf_class One = 1.0;
    mpf_class ctemp, stemp, temp;
    mpackint info;
    mpackint i, j;

    info = 0;
    if (!(Mlsame_gmp(side, "L") || Mlsame_gmp(side, "R")))
	info = 1;
    else if (!(Mlsame_gmp(pivot, "V") || Mlsame_gmp(pivot, "T")
	    || Mlsame_gmp(pivot, "B")))
	info = 2;
    else if (!(Mlsame_gmp(direct, "F") || Mlsame_gmp(direct, "B")))
	info = 3;
    else if (m < 0)
	info = 4;
    else if (n < 0)
	info = 5;
    else if (lda < max((mpackint)1, m))
	info = 9;
    if (info != 0) {
	Mxerbla_gmp("Rlasr ", info);
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0) {
	return;
    }

    if (Mlsame_gmp(side, "L")) {
//Form  P * A
	if (Mlsame_gmp(pivot, "V")) {
	    if (Mlsame_gmp(direct, "F")) {
		for (j = 0; j < m - 1; j++) {
		    ctemp = c[j];
		    stemp = s[j];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < n; i++) {
			    temp = A[(j + 1) + i * lda];
			    A[(j + 1) + i * lda] = ctemp * temp - stemp *
				A[j + i * lda];
			    A[j + i * lda] =
				stemp * temp + ctemp * A[j + i * lda];
			}
		    }
		}
	    } else if (Mlsame_gmp(direct, "B")) {
		for (j = m - 2; j >= 0; j--) {
		    ctemp = c[j];
		    stemp = s[j];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < n; i++) {
			    temp = A[(j + 1) + i * lda];
			    A[(j + 1) + i * lda] = ctemp * temp - stemp *
				A[j + i * lda];
			    A[j + i * lda] = stemp * temp + ctemp * A[j
				+ i * lda];
			}
		    }
		}
	    }
	}

	else if (Mlsame_gmp(pivot, "T")) {
	    if (Mlsame_gmp(direct, "F")) {
		for (j = 1; j < m; j++) {
		    ctemp = c[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < n; i++) {
			    temp = A[j + i * lda];
			    A[j + i * lda] = ctemp * temp - stemp * A[i * lda];
			    A[i * lda] = stemp * temp + ctemp * A[i * lda];
			}
		    }
		}
	    } else if (Mlsame_gmp(direct, "B")) {
		for (j = m - 1; j >= 1; j--) {
		    ctemp = c[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < n; i++) {
			    temp = A[j + i * lda];
			    A[j + i * lda] = ctemp * temp - stemp * A[i * lda];
			    A[i * lda] = stemp * temp + ctemp * A[i * lda];
			}
		    }
		}
	    }
	}

	else if (Mlsame_gmp(pivot, "B")) {
	    if (Mlsame_gmp(direct, "F")) {
		for (j = 0; j < m - 1; j++) {
		    ctemp = c[j];
		    stemp = s[j];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < n; i++) {
			    temp = A[j + i * lda];
			    A[j + i * lda] = stemp * A[(m - 1) + i * lda]
				+ ctemp * temp;
			    A[(m - 1) + i * lda] =
				ctemp * A[(m - 1) + i * lda] - stemp * temp;
			}
		    }
		}
	    } else if (Mlsame_gmp(direct, "B")) {
		for (j = m - 2; j >= 0; j--) {
		    ctemp = c[j];
		    stemp = s[j];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < n; i++) {
			    temp = A[j + i * lda];
			    A[j + i * lda] = stemp * A[(m - 1) + i * lda]
				+ ctemp * temp;
			    A[(m - 1) + i * lda] =
				ctemp * A[(m - 1) + i * lda] - stemp * temp;
			}
		    }
		}
	    }
	}
    }

    else if (Mlsame_gmp(side, "R")) {
//Form A * P'
	if (Mlsame_gmp(pivot, "V")) {
	    if (Mlsame_gmp(direct, "F")) {
		for (j = 0; j < n - 1; j++) {
		    ctemp = c[j];
		    stemp = s[j];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < m; i++) {
			    temp = A[i + (j + 1) * lda];
			    A[i + (j + 1) * lda] =
				ctemp * temp - stemp * A[i + j * lda];
			    A[i + j * lda] =
				stemp * temp + ctemp * A[i + j * lda];
			}
		    }
		}
	    } else if (Mlsame_gmp(direct, "B")) {
		for (j = n - 2; j >= 0; j--) {
		    ctemp = c[j];
		    stemp = s[j];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < m; i++) {
			    temp = A[i + (j + 1) * lda];
			    A[i + (j + 1) * lda] =
				ctemp * temp - stemp * A[i + j * lda];
			    A[i + j * lda] =
				stemp * temp + ctemp * A[i + j * lda];
			}
		    }
		}
	    }
	} else if (Mlsame_gmp(pivot, "T")) {
	    if (Mlsame_gmp(direct, "F")) {
		for (j = 1; j < n; j++) {
		    ctemp = c[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < m; i++) {
			    temp = A[i + j * lda];
			    A[i + j * lda] = ctemp * temp - stemp * A[i];
			    A[i] = stemp * temp + ctemp * A[i];
			}
		    }
		}
	    } else if (Mlsame_gmp(direct, "B")) {
		for (j = n - 1; j >= 1; j--) {
		    ctemp = c[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < m; i++) {
			    temp = A[i + j * lda];
			    A[i + j * lda] = ctemp * temp - stemp * A[i];
			    A[i] = stemp * temp + ctemp * A[i];
			}
		    }
		}
	    }
	} else if (Mlsame_gmp(pivot, "B")) {
	    if (Mlsame_gmp(direct, "F")) {
		for (j = 0; j < n - 1; j++) {
		    ctemp = c[j];
		    stemp = s[j];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < m; i++) {
			    temp = A[i + j * lda];
			    A[i + j * lda] = stemp * A[i + (n - 1) * lda]
				+ ctemp * temp;
			    A[i + (n - 1) * lda] =
				ctemp * A[i + (n - 1) * lda] - stemp * temp;
			}
		    }
		}
	    } else if (Mlsame_gmp(direct, "B")) {
		for (j = n - 2; j >= 0; j--) {
		    ctemp = c[j];
		    stemp = s[j];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < m; i++) {
			    temp = A[i + j * lda];
			    A[i + j * lda] = stemp * A[i + (n - 1) * lda]
				+ ctemp * temp;
			    A[i + (n - 1) * lda] =
				ctemp * A[i + (n - 1) * lda] - stemp * temp;
			}
		    }
		}
	    }
	}
    }
    return;
}
