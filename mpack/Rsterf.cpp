/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rsterf.cpp,v 1.8 2009/09/26 02:21:32 nakatamaho Exp $ 
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
#include <stdio.h> //for untested part
void
Rsterf(mpackint n, mpf_class * d, mpf_class * e, mpackint *info)
{
    mpf_class Zero = 0.0, One = 1.0, Two = 2.0, Three = 3.0;
    mpf_class sigma;
    mpf_class eps, eps2;
    mpf_class safmin, safmax, ssfmax, ssfmin, anorm;
    mpf_class rte, rt1, rt2, s, c, r, oldc, oldgam, gamma, p, bb, alpha;

    mpackint nmaxit;
    mpackint iscale;
    mpackint l1, jtot;
    mpackint i, l, m;
    mpackint lsv, lend, lendsv;

    *info = 0;
//Quick return if possible
    if (n < 0) {
	*info = -1;
	Mxerbla_gmp("Rsterf", -(*info));
	return;
    }
    if (n <= 1)
	return;
//Determine the unit roundoff for this environment.
    eps = Rlamch_gmp("E");
    eps2 = eps * eps;
    safmin = Rlamch_gmp("S");
    safmax = One / safmin;
    ssfmax = sqrt(safmax) / Three;
    ssfmin = sqrt(safmin) / eps2;
//Compute the eigenvalues of the tridiagonal matrix.
    nmaxit = n * 30;
    sigma = Zero;
    jtot = 0;
//Determine where the matrix splits and choose QL or QR iteration
//for each block, according to whether top or bottom diagonal
//element is smaller.
    l1 = 1;
  L10:
    if (l1 > n) {
	goto L170;
    }
    if (l1 > 1) {
	e[l1 - 2] = Zero;
    }
    for (m = l1; m <= n - 1; m++) {
	if (abs(e[m - 1]) <= sqrt(abs(d[m - 1])) * sqrt(abs(d[m])) * eps) {
	    e[m - 1] = Zero;
	    goto L30;
	}
    }
    m = n;
  L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
	goto L10;
    }
//Scale submatrix in rows and columns L to LEND
    anorm = Rlanst("I", lend - l + 1, &d[l - 1], &e[l - 1]);
    iscale = 0;
    if (anorm > ssfmax) {
	printf("XXX not tested #1\n");
	iscale = 1;
	Rlascl("G", 0, 0, anorm, ssfmax, lend - l + 1, 1, &d[l - 1], n, info);
	Rlascl("G", 0, 0, anorm, ssfmax, lend - l, 1, &e[l - 1], n, info);
    } else if (anorm < ssfmin) {
	printf("XXX not tested #2\n");
	iscale = 2;
	Rlascl("G", 0, 0, anorm, ssfmin, lend - l + 1, 1, &d[l - 1], n, info);
	Rlascl("G", 0, 0, anorm, ssfmin, lend - l, 1, &e[l - 1], n, info);
    }
    for (i = l; i <= lend - 1; i++) {
	e[i - 1] = e[i - 1] * e[i - 1];
    }
//Choose between QL and QR iteration
    if (abs(d[lend - 1]) < abs(d[l - 1])) {
	lend = lsv;
	l = lendsv;
    }
    if (lend >= l) {
//QL Iteration
//Look for small subdiagonal element.
      L50:
	if (l != lend) {
	    for (m = l; m <= lend - 1; m++) {
		if (abs(e[m - 1]) <= eps2 * abs(d[m - 1] * d[m])) {
		    goto L70;
		}
	    }
	}
	m = lend;
      L70:
	if (m < lend) {
	    e[m - 1] = Zero;
	}
	p = d[l - 1];
	if (m == l) {
	    goto L90;
	}
//If remaining matrix is 2 by 2, use DLAE2 to compute its
//eigenvalues.
	if (m == l + 1) {
	    rte = sqrt(e[l - 1]);
	    Rlae2(d[l - 1], rte, d[l], &rt1, &rt2);
	    d[l - 1] = rt1;
	    d[l] = rt2;
	    e[l - 1] = Zero;
	    l = l + 2;
	    if (l <= lend) {
		goto L50;
	    }
	    goto L150;
	}
	if (jtot == nmaxit) {
	    goto L150;
	}
	jtot++;
//Form shift.
	rte = sqrt(e[l - 1]);
	sigma = (d[l] - p) / (rte * Two);
	r = Rlapy2(sigma, One);
	sigma = p - rte / (sigma + Msign(r, sigma));
	c = One;
	s = Zero;
	gamma = d[m - 1] - sigma;
	p = gamma * gamma;
//Inner loop 
	for (i = m - 1; i >= l; i--) {
	    bb = e[i - 1];
	    r = p + bb;
	    if (i != m - 1) {
		e[i] = s * r;
	    }
	    oldc = c;
	    c = p / r;
	    s = bb / r;
	    oldgam = gamma;
	    alpha = d[i - 1];
	    gamma = c * (alpha - sigma) - s * oldgam;
	    d[i] = oldgam + (alpha - gamma);
	    if (c != Zero) {
		p = gamma * gamma / c;
	    } else {
		p = oldc * bb;
	    }
	}
	e[l - 1] = s * p;
	d[l - 1] = sigma + gamma;
	goto L50;
//Eigenvalue found.
      L90:
	d[l - 1] = p;
	l++;
	if (l <= lend) {
	    goto L50;
	}
	goto L150;
    } else {
//QR Iteration
//Look for small superdiagonal element.
      L100:
	for (m = l; m >= lend + 1; m--) {
	    if (abs(e[m - 2]) <= eps2 * abs(d[m - 1] * d[m - 2])) {
		goto L120;
	    }
	}
	m = lend;
      L120:
	if (m > lend) {
	    e[m - 2] = Zero;
	}
	p = d[l - 1];
	if (m == l) {
	    goto L140;
	}
//If remaining matrix is 2 by 2, use DLAE2 to compute its
//eigenvalues.
	if (m == l - 1) {
	    rte = sqrt(e[l - 2]);
	    Rlae2(d[l - 1], rte, d[l - 2], &rt1, &rt2);
	    d[l - 1] = rt1;
	    d[l - 2] = rt2;
	    e[l - 2] = Zero;
	    l = l - 2;
	    if (l >= lend) {
		goto L100;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	jtot++;
//Form shift.
	rte = sqrt(e[l - 2]);
	sigma = (d[l - 2] - p) / (rte * Two);
	r = Rlapy2(sigma, One);
	sigma = p - rte / (sigma + Msign(r, sigma));

	c = One;
	s = Zero;
	gamma = d[m - 1] - sigma;
	p = gamma * gamma;
//Inner loop
	for (i = m; i <= l - 1; i++) {
	    bb = e[i - 1];
	    r = p + bb;
	    if (i != m) {
		e[i - 2] = s * r;
	    }
	    oldc = c;
	    c = p / r;
	    s = bb / r;
	    oldgam = gamma;
	    alpha = d[i];
	    gamma = c * (alpha - sigma) - s * oldgam;
	    d[i - 1] = oldgam + (alpha - gamma);
	    if (c != Zero) {
		p = gamma * gamma / c;
	    } else {
		p = oldc * bb;
	    }
	}

	e[l - 2] = s * p;
	d[l - 1] = sigma + gamma;
	goto L100;
//Eigenvalue found.
      L140:
	d[l - 1] = p;

	l--;
	if (l >= lend) {
	    goto L100;
	}
	goto L150;
    }
//Undo scaling if necessary
  L150:
    if (iscale == 1) {
	Rlascl("G", 0, 0, ssfmax, anorm, lendsv - lsv + 1, 1, &d[lsv - 1], n,
	    info);
    }
    if (iscale == 2) {
	Rlascl("G", 0, 0, ssfmin, anorm, lendsv - lsv + 1, 1, &d[lsv - 1], n,
	    info);
    }
//Check for no convergence to an eigenvalue after a total
//of N*MAXIT iterations.
    if (jtot < nmaxit) {
	goto L10;
    }
    for (i = 1; i <= n - 1; i++) {
	if (e[i - 1] != Zero) {
	    ++(*info);
	}
    }
    return;

//Sort eigenvalues in increasing order.
  L170:
    Rlasrt("I", n, &d[0], info);
    return;
}
