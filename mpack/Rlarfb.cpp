/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 *
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rlarfb.cpp,v 1.2 2009/09/12 07:59:10 nakatamaho Exp $ 
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
 * Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.
 * 
 * $COPYRIGHT$
 * 
 * Additional copyrights may follow
 * 
 * $HEADER$
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * - Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer listed in this
 * license in the documentation and/or other materials provided with the
 * distribution.
 * 
 * - Neither the name of the copyright holders nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <mblas_gmp.h>
#include <mlapack_gmp.h>

void
Rlarfb(const char *side, const char *trans, const char *direct,
    const char *storev, mpackint m, mpackint n, mpackint k, mpf_class * V, mpackint ldv,
    mpf_class * T, mpackint ldt, mpf_class * C, mpackint ldc, mpf_class * work,
    mpackint ldwork)
{
    mpackint i, j;
    mpf_class One = 1.0;
    mpf_class mOne = -1.0;
    char transt;

    //Quick return if possible
    if (m <= 0 || n <= 0)
	return;

    if (Mlsame_gmp(trans, "N")) {
	transt = 'T';
    } else {
	transt = 'N';
    }

    if (Mlsame_gmp(storev, "C")) {
	if (Mlsame_gmp(direct, "F")) {

//Let V = (V1) (first K rows)
//        (V2)
// where V1 is unit lower triangular.
	    if (Mlsame_gmp(side, "L")) {

//Form H * C or H ' * C  where  C = ( C1 )
//                                  ( C2 )
// W: = C ' * V  =  (C1' * V1 + C2 '*V2)  (stored in WORK)
// W: = C1 '
		for (j = 0; j < k; j++) {
		    Rcopy(n, &C[j], ldc, &work[j * ldwork], 1);
		}
//W: = W * V1
		Rtrmm("Right", "Lower", "No transpose", "Unit", n, k, One,
		    &V[0], ldv, &work[0], ldwork);
		if (m > k) {
//W: = W + C2 '*V2
		    Rgemm("Transpose", "No transpose", n, k, m - k, One,
			&C[k], ldc, &V[k], ldv, One, &work[0], ldwork);
		}
//W: = W * T '  or  W * T
		Rtrmm("Right", "Upper", &transt, "Non-unit", n, k, One, &T[0],
		    ldt, &work[0], ldwork);
//C: = C - V * W '
		if (m > k) {
//C2: = C2 - V2 * W '
		    Rgemm("No transpose", "Transpose", m - k, n, k, mOne,
			&V[k], ldv, &work[0], ldwork, One, &C[k], ldc);
		}
//W: = W * V1 '
		Rtrmm("Right", "Lower", "Transpose", "Unit", n, k, One, &V[0],
		    ldv, &work[0], ldwork);

//C1: = C1 - W '
		for (j = 0; j < k; j++) {
		    for (i = 0; i < n; i++) {
			C[j + i * ldc] -= work[i + j * ldwork];
		    }
		}
	    } else if (Mlsame_gmp(side, "R")) {
//Form C * H or C * H '  where  C = ( C1  C2 )
//W: = C * V = (C1 * V1 + C2 * V2) (stored in WORK)
//W: = C1
		for (j = 0; j < k; j++) {
		    Rcopy(m, &C[j * ldc], 1, &work[j * ldwork], 1);
		}
//W: = W * V1
		Rtrmm("Right", "Lower", "No transpose", "Unit", m, k, One,
		    &V[0], ldv, &work[0], ldwork);
		if (n > k) {
//W: = W + C2 * V2
		    Rgemm("No transpose", "No transpose", m, k, n - k, One,
			&C[k * ldc], ldc, &V[k], ldv, One, &work[0], ldwork);
		}
//W: = W * T or W * T '
		Rtrmm("Right", "Upper", trans, "Non-unit", m, k, One, &T[0],
		    ldt, &work[0], ldwork);
//C: = C - W * V '
		if (n > k) {
//C2: = C2 - W * V2'
		    Rgemm("No transpose", "Transpose", m, n - k, k, mOne,
			&work[0], ldwork, &V[k], ldv, One, &C[k * ldc], ldc);
		}
//W: = W * V1 '
		Rtrmm("Right", "Lower", "Transpose", "Unit", m, k, One, &V[0],
		    ldv, &work[0], ldwork);
//C1: = C1 - W
		for (j = 0; j < k; j++) {
		    for (i = 0; i < m; i++) {
			C[i + j * ldc] -= work[i + j * ldwork];
		    }
		}
	    }
	} else {
//Let V = (V1)
//        (V2) (last K rows)
// where V2 is unit upper triangular.
	    if (Mlsame_gmp(side, "L")) {
//Form H * C or H ' * C  where  C = ( C1 )
//                                  ( C2 )
//W: = C ' * V  =  (C1' * V1 + C2 '*V2)  (stored in WORK)
//W: = C2 '
		for (j = 0; j < k; j++) {
		    Rcopy(n, &C[m - k + j], ldc, &work[j * ldwork], 1);
		}
//W: = W * V2
		Rtrmm("Right", "Upper", "No transpose", "Unit", n, k, One,
		    &V[m - k], ldv, &work[0], ldwork);
		if (m > k) {
//W: = W + C1 '*V1
		    Rgemm("Transpose", "No transpose", n, k, m - k, One,
			&C[0], ldc, &V[0], ldv, One, &work[0], ldwork);
		}
//W: = W * T '  or  W * T
		Rtrmm("Right", "Lower", &transt, "Non-unit", n, k, One, &T[0],
		    ldt, &work[0], ldwork);
//C: = C - V * W '
		if (m > k) {
//C1:= C1 - V1 * W '
		    Rgemm("No transpose", "Transpose", m - k, n, k, mOne,
			&V[0], ldv, &work[0], ldwork, One, &C[0], ldc);
		}
//W: = W * V2 '
		Rtrmm("Right", "Upper", "Transpose", "Unit", n, k, One,
		    &V[m - k], ldv, &work[0], ldwork);
//C2:= C2 - W '
		for (j = 0; j < k; j++) {
		    for (i = 0; i < n; i++) {
			C[m - k + j + i * ldc] -= work[i + j * ldwork];
		    }
		}
	    } else if (Mlsame_gmp(side, "R")) {
//Form C * H or C * H '  where  C = ( C1  C2 )
// W: = C * V = (C1 * V1 + C2 * V2) (stored in WORK)
// W: = C2
		for (j = 0; j < k; j++) {
		    Rcopy(m, &C[(n - k + j) * ldc], 1, &work[j * ldwork], 1);
		}
//W:= W * V2
		Rtrmm("Right", "Upper", "No transpose", "Unit", m, k, One,
		    &V[n - k], ldv, &work[0], ldwork);
		if (n > k) {
//W:= W + C1 * V1
		    Rgemm("No transpose", "No transpose", m, k, n - k, One,
			&C[0], ldc, &V[0], ldv, One, &work[0], ldwork);
		}
//W:= W * T or W * T
		Rtrmm("Right", "Lower", trans, "Non-unit", m, k, One, &T[0],
		    ldt, &work[0], ldwork);
//C:= C - W * V '
		if (n > k) {
//C1:= C1 - W * V1 '
		    Rgemm("No transpose", "Transpose", m, n - k, k, mOne,
			&work[0], ldwork, &V[0], ldv, One, &C[0], ldc);
		}
//W: = W * V2 '
		Rtrmm("Right", "Upper", "Transpose", "Unit", m, k, One,
		    &V[n - k], ldv, &work[0], ldwork);
//C2:= C2 - W
		for (j = 0; j < k; j++) {
		    for (i = 0; i < m; i++) {
			C[i + (n - k + j) * ldc] -= work[i + j * ldwork];
		    }
		}
	    }
	}
    } else if (Mlsame_gmp(storev, "R")) {
	if (Mlsame_gmp(direct, "F")) {
//Let V = (V1 V2) (V1:first K columns)
//where V1 is unit upper triangular.

	    if (Mlsame_gmp(side, "L")) {
//Form H * C or H ' * C  where  C = ( C1 )
//                                  ( C2 )
// W:= C ' * V' = (C1 '*V1' + C2 '*V2') (stored in WORK)
// W:= C1 '
		for (j = 0; j < k; j++) {
		    Rcopy(n, &C[j], ldc, &work[j * ldwork], 1);
		}
//W:= W * V1 '
		Rtrmm("Right", "Upper", "Transpose", "Unit", n, k, One, &V[0],
		    ldv, &work[0], ldwork);
		if (m > k) {
//W:= W + C2 '*V2'
		    Rgemm("Transpose", "Transpose", n, k, m - k, One,
			&C[k], ldc, &V[k * ldv], ldv, One, &work[0], ldwork);
		}
//W:= W * T '  or  W * T
		Rtrmm("Right", "Upper", &transt, "Non-unit", n, k, One,
		    &T[0], ldt, &work[0], ldwork);
//C:= C - V ' * W'
		if (m > k) {
//C2:= C2 - V2 ' * W'
		    Rgemm("Transpose", "Transpose", m - k, n, k, mOne,
			&V[k * ldv], ldv, &work[0], ldwork, One, &C[k], ldc);
		}
//W:= W * V1
		Rtrmm("Right", "Upper", "No transpose", "Unit", n, k, One,
		    &V[0], ldv, &work[0], ldwork);
//C1:= C1 - W '
		for (j = 0; j < k; j++) {
		    for (i = 0; i < n; i++) {
			C[j + i * ldc] -= work[i + j * ldwork];
		    }
		}
	    } else if (Mlsame_gmp(side, "R")) {
//Form C * H or C * H '  where  C = ( C1  C2 )
// W:= C * V '  =  (C1*V1' + C2 * V2 ')  (stored in WORK)
// W:= C1
		for (j = 0; j < k; j++) {
		    Rcopy(m, &C[j * ldc], 1, &work[j * ldwork], 1);
		}
//W:= W * V1 '
		Rtrmm("Right", "Upper", "Transpose", "Unit", m, k, One, &V[0],
		    ldv, &work[0], ldwork);
		if (n > k) {
//W:= W + C2 * V2 '
		    Rgemm("No transpose", "Transpose", m, k, n - k, One,
			&C[k * ldc], ldc, &V[k * ldv],
			ldv, One, &work[0], ldwork);
		}
//W:= W * T or W * T '
		Rtrmm("Right", "Upper", trans, "Non-unit", m, k, One, &T[0],
		    ldt, &work[0], ldwork);
//C:= C - W * V
		if (n > k) {
//C2:= C2 - W * V2
		    Rgemm("No transpose", "No transpose", m, n - k, k, mOne,
			&work[0], ldwork, &V[k * ldv], ldv, One,
			&C[k * ldc], ldc);
		}
//W:= W * V1
		Rtrmm("Right", "Upper", "No transpose", "Unit", m, k, One,
		    &V[0], ldv, &work[0], ldwork);
//C1:= C1 - W
		for (j = 0; j < k; j++) {
		    for (i = 0; i < m; i++) {
			C[i + j * ldc] -= work[i + j * ldwork];
		    }
		}
	    }
	} else {
//Let V = (V1 V2) (V2:last K columns)
// where V2 is unit lower triangular.
	    if (Mlsame_gmp(side, "L")) {
//Form H * C or H ' * C  where  C = ( C1 )
//                                  ( C2 )
//W:= C ' * V' = (C1 '*V1' + C2 '*V2') (stored in WORK)
//W:= C2 '
		for (j = 0; j < k; j++) {
		    Rcopy(n, &C[m - k + j], ldc, &work[j * ldwork], 1);
		}

//W:= W * V2 '
		Rtrmm("Right", "Lower", "Transpose", "Unit", n, k, One,
		    &V[(m - k) * ldv], ldv, &work[0], ldwork);

		if (m > k) {

//W:= W + C1 '*V1'
		    Rgemm("Transpose", "Transpose", n, k, m - k, One, &C[0],
			ldc, &V[0], ldv, One, &work[0], ldwork);
		}
//W:= W * T '  or  W * T
		Rtrmm("Right", "Lower", &transt, "Non-unit", n, k, One, &T[0],
		    ldt, &work[0], ldwork);
//C:= C - V ' * W'
		if (m > k) {

//C1:= C1 - V1 ' * W'
		    Rgemm("Transpose", "Transpose", m - k, n, k, mOne, &V[0],
			ldv, &work[0], ldwork, One, &C[0], ldc);
		}
//W:= W * V2
		Rtrmm("Right", "Lower", "No transpose", "Unit", n, k, One,
		    &V[(m - k) * ldv], ldv, &work[0], ldwork);
//C2:= C2 - W '
		for (j = 0; j < k; j++) {
		    for (i = 0; i < n; i++) {
			C[m - k + j + i * ldc] -= work[i + j * ldwork];
		    }
		}
	    } else if (Mlsame_gmp(side, "R")) {
//Form C * H or C * H '  where  C = ( C1  C2 )
// W:= C * V '  =  (C1*V1' + C2 * V2 ')  (stored in WORK)
// W:= C2
		for (j = 0; j < k; j++) {
		    Rcopy(m, &C[(n - k + j) * ldc], 1, &work[j * ldwork], 1);
		}
//W: = W * V2 '
		Rtrmm("Right", "Lower", "Transpose", "Unit", m, k, One,
		    &V[(n - k) * ldv], ldv, &work[0], ldwork);
		if (n > k) {
//W:= W + C1 * V1 '
		    Rgemm("No transpose", "Transpose", m, k, n - k, One, &C[0],
			ldc, &V[0], ldv, One, &work[0], ldwork);
		}
//W:= W * T or W * T '
		Rtrmm("Right", "Lower", trans, "Non-unit", m, k, One, &T[0],
		    ldt, &work[0], ldwork);
//C:= C - W * V
		if (n > k) {
//C1:= C1 - W * V1
		    Rgemm("No transpose", "No transpose", m, n - k, k, mOne,
			&work[0], ldwork, &V[0], ldv, One, &C[0], ldc);
		}
//W:=W * V2
		Rtrmm("Right", "Lower", "No transpose", "Unit", m, k, One,
		    &V[(n - k) * ldv], ldv, &work[0], ldwork);
//C1: = C1 - W
		for (j = 0; j < k; j++) {
		    for (i = 0; i < m; i++) {
			C[i + (n - k + j) * ldc] -= work[i + j * ldwork];
		    }
		}
	    }
	}
    }
    return;
}
