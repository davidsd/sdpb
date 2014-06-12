/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: mlapack_gmp.h,v 1.6 2009/09/22 20:27:18 nakatamaho Exp $ 
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

#ifndef _MLAPACK_GMP_H_
#define _MLAPACK_GMP_H_

/* this is a subset of mpack for SDPA-GMP only */
/* http://mplapack.sourceforge.net/ */

/* mlapack prototypes */
void Rsteqr(const char *compz, mpackint n, mpf_class * d, mpf_class * e,
    mpf_class * Z, mpackint ldz, mpf_class * work, mpackint *info);
void
    Rsyev(const char *jobz, const char *uplo, mpackint n, mpf_class * A,
    mpackint lda, mpf_class * w, mpf_class * work, mpackint *lwork, mpackint *info);
void Rpotrf(const char *uplo, mpackint n, mpf_class * A, mpackint lda, mpackint *info);
mpackint iMlaenv_gmp(mpackint ispec, const char *name, const char *opts, mpackint n1, mpackint n2,
    mpackint n3, mpackint n4);
mpf_class Rlamch_gmp(const char *cmach);
mpf_class Rlansy(const char *norm, const char *uplo, mpackint n, mpf_class * A,
    mpackint lda, mpf_class * work);
void Rlascl(const char *type, mpackint kl, mpackint ku, mpf_class cfrom, mpf_class cto,
    mpackint m, mpackint n, mpf_class * A, mpackint lda, mpackint *info);
void Rsytrd(const char *uplo, mpackint n, mpf_class * A, mpackint lda, mpf_class * d,
    mpf_class * e, mpf_class * tau, mpf_class * work, mpackint lwork, mpackint *info);
void Rsytd2(const char *uplo, mpackint n, mpf_class * A, mpackint lda, mpf_class * d,
    mpf_class * e, mpf_class * tau, mpackint *info);
mpf_class Rlanst(const char *norm, mpackint n, mpf_class * d, mpf_class * e);
void Rlae2(mpf_class a, mpf_class b, mpf_class c, mpf_class * rt1,
    mpf_class * rt2);
mpf_class Rlapy2(mpf_class x, mpf_class y);
void Rlasrt(const char *id, mpackint n, mpf_class * d, mpackint *info);
void Rorgql(mpackint m, mpackint n, mpackint k, mpf_class * A, mpackint lda, mpf_class * tau,
    mpf_class * work, mpackint lwork, mpackint *info);
void Rorgqr(mpackint m, mpackint n, mpackint k, mpf_class * A, mpackint lda, mpf_class * tau,
    mpf_class * work, mpackint lwork, mpackint *info);
void Rlarfg(mpackint N, mpf_class * alpha, mpf_class * x, mpackint incx,
    mpf_class * tau);
void Rlassq(mpackint n, mpf_class * x, mpackint incx, mpf_class * scale,
    mpf_class * sumsq);
void Rorg2l(mpackint m, mpackint n, mpackint k, mpf_class * A, mpackint lda, mpf_class * tau,
    mpf_class * work, mpackint *info);
void Rlarft(const char *direct, const char *storev, mpackint n, mpackint k,
    mpf_class * v, mpackint ldv, mpf_class * tau, mpf_class * t, mpackint ldt);
void Rlarfb(const char *side, const char *trans, const char *direct,
    const char *storev, mpackint m, mpackint n, mpackint k, mpf_class * V, mpackint ldv,
    mpf_class * T, mpackint ldt, mpf_class * C, mpackint ldc, mpf_class * work,
    mpackint ldwork);
void Rorg2r(mpackint m, mpackint n, mpackint k, mpf_class * A, mpackint lda, mpf_class * tau,
    mpf_class * work, mpackint *info);
void Rlarf(const char *side, mpackint m, mpackint n, mpf_class * v, mpackint incv,
    mpf_class tau, mpf_class * C, mpackint ldc, mpf_class * work);
void Rpotf2(const char *uplo, mpackint n, mpf_class * A, mpackint lda, mpackint *info);
void Rlaset(const char *uplo, mpackint m, mpackint n, mpf_class alpha, mpf_class beta,
    mpf_class * A, mpackint lda);
void Rlaev2(mpf_class a, mpf_class b, mpf_class c, mpf_class * rt1,
    mpf_class * rt2, mpf_class * cs1, mpf_class * sn1);
void Rlasr(const char *side, const char *pivot, const char *direct, mpackint m,
    mpackint n, mpf_class * c, mpf_class * s, mpf_class * A, mpackint lda);
void Rlartg(mpf_class f, mpf_class g, mpf_class * cs, mpf_class * sn,
    mpf_class * r);
void Rlatrd(const char *uplo, mpackint n, mpackint nb, mpf_class * A, mpackint lda, mpf_class * e, mpf_class * tau, mpf_class * w, mpackint ldw);
void Rsterf(mpackint n, mpf_class * d, mpf_class * e, mpackint *info);
void Rorgtr(const char *uplo, mpackint n, mpf_class * a, mpackint lda, mpf_class * tau,
    mpf_class * work, mpackint lwork, mpackint *info);
#endif
