/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rlarfg.cpp,v 1.4 2009/09/26 02:21:32 nakatamaho Exp $ 
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
#include <stdio.h> //for debugging
void
Rlarfg(mpackint N, mpf_class * alpha, mpf_class * x, mpackint incx, mpf_class * tau)
{
    mpf_class xnorm;
    mpf_class Zero = 0.0, One = 1.0;
    mpf_class beta;
    mpf_class safmin;
    mpf_class rsafmn;
    mpackint knt;

    if (N <= 1) {
	*tau = 0.0;
	return;
    }
    xnorm = Rnrm2(N - 1, x, incx);
//H  =  I
    if (xnorm == 0.0) {
	*tau = 0.0;
    } else {
	beta = -1.0 * Msign(Rlapy2(*alpha, xnorm), *alpha);
	safmin = Rlamch_gmp("S") / Rlamch_gmp("E");

//XNORM, BETA may be inaccurate; scale X and recompute them
	if (abs(beta) < safmin) {
	    fprintf(stderr, "# Rlarfg: 1: XXX not very well tested\n");
	    rsafmn = One / safmin;
	    knt = 0;
	    while (abs(beta) < safmin) {
		knt++;
		Rscal(N - 1, rsafmn, x, incx);
		beta = beta * rsafmn;
		*alpha = *alpha * rsafmn;
	    }

//New BETA is at most 1, at least SAFMIN
	    xnorm = Rnrm2(N - 1, x, incx);
	    beta = -1.0 * Msign(Rlapy2(*alpha, xnorm), *alpha);
	    *tau = (beta - *alpha) / beta;
	    Rscal(N - 1, One / (*alpha - beta), x, incx);

//If ALPHA is subnormal, it may lose relative accuracy
	    *alpha = beta;
	    for (mpackint j = 0; j < knt; j++) {
		*alpha = *alpha * safmin;
	    }
	} else {
	    *tau = (beta - *alpha) / beta;
	    Rscal(N - 1, One / (*alpha - beta), x, incx);
	    *alpha = beta;
	}
    }
}
