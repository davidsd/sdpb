/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: mlapack_mpfr.h,v 1.6 2009/09/22 20:27:18 nakatamaho Exp $ 
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

#ifndef _MLAPACK_MPFR_H_
#define _MLAPACK_MPFR_H_

/* this is a subset of mpack for SDPA-GMP only */
/* http://mplapack.sourceforge.net/ */

/* mlapack prototypes */
void Rsteqr(const char *compz, mpackint n, mpreal * d, mpreal * e,
    mpreal * Z, mpackint ldz, mpreal * work, mpackint *info);
void
    Rsyev(const char *jobz, const char *uplo, mpackint n, mpreal * A,
    mpackint lda, mpreal * w, mpreal * work, mpackint lwork, mpackint *info);
void Rpotrf(const char *uplo, mpackint n, mpreal * A, mpackint lda, mpackint *info);
mpackint iMlaenv_mpfr(mpackint ispec, const char *name, const char *opts, mpackint n1, mpackint n2,
    mpackint n3, mpackint n4);
mpreal Rlamch_mpfr(const char *cmach);
mpreal Rlansy(const char *norm, const char *uplo, mpackint n, mpreal * A,
    mpackint lda, mpreal * work);
void Rlascl(const char *type, mpackint kl, mpackint ku, mpreal cfrom, mpreal cto,
    mpackint m, mpackint n, mpreal * A, mpackint lda, mpackint *info);
void Rsytrd(const char *uplo, mpackint n, mpreal * A, mpackint lda, mpreal * d,
    mpreal * e, mpreal * tau, mpreal * work, mpackint lwork, mpackint *info);
void Rsytd2(const char *uplo, mpackint n, mpreal * A, mpackint lda, mpreal * d,
    mpreal * e, mpreal * tau, mpackint *info);
mpreal Rlanst(const char *norm, mpackint n, mpreal * d, mpreal * e);
void Rlae2(mpreal a, mpreal b, mpreal c, mpreal * rt1,
    mpreal * rt2);
mpreal Rlapy2(mpreal x, mpreal y);
void Rlasrt(const char *id, mpackint n, mpreal * d, mpackint *info);
void Rorgql(mpackint m, mpackint n, mpackint k, mpreal * A, mpackint lda, mpreal * tau,
    mpreal * work, mpackint lwork, mpackint *info);
void Rorgqr(mpackint m, mpackint n, mpackint k, mpreal * A, mpackint lda, mpreal * tau,
    mpreal * work, mpackint lwork, mpackint *info);
void Rlarfg(mpackint N, mpreal * alpha, mpreal * x, mpackint incx,
    mpreal * tau);
void Rlassq(mpackint n, mpreal * x, mpackint incx, mpreal * scale,
    mpreal * sumsq);
void Rorg2l(mpackint m, mpackint n, mpackint k, mpreal * A, mpackint lda, mpreal * tau,
    mpreal * work, mpackint *info);
void Rlarft(const char *direct, const char *storev, mpackint n, mpackint k,
    mpreal * v, mpackint ldv, mpreal * tau, mpreal * t, mpackint ldt);
void Rlarfb(const char *side, const char *trans, const char *direct,
    const char *storev, mpackint m, mpackint n, mpackint k, mpreal * V, mpackint ldv,
    mpreal * T, mpackint ldt, mpreal * C, mpackint ldc, mpreal * work,
    mpackint ldwork);
void Rorg2r(mpackint m, mpackint n, mpackint k, mpreal * A, mpackint lda, mpreal * tau,
    mpreal * work, mpackint *info);
void Rlarf(const char *side, mpackint m, mpackint n, mpreal * v, mpackint incv,
    mpreal tau, mpreal * C, mpackint ldc, mpreal * work);
void Rpotf2(const char *uplo, mpackint n, mpreal * A, mpackint lda, mpackint *info);
void Rlaset(const char *uplo, mpackint m, mpackint n, mpreal alpha, mpreal beta,
    mpreal * A, mpackint lda);
void Rlaev2(mpreal a, mpreal b, mpreal c, mpreal * rt1,
    mpreal * rt2, mpreal * cs1, mpreal * sn1);
void Rlasr(const char *side, const char *pivot, const char *direct, mpackint m,
    mpackint n, mpreal * c, mpreal * s, mpreal * A, mpackint lda);
void Rlartg(mpreal f, mpreal g, mpreal * cs, mpreal * sn,
    mpreal * r);
void Rlatrd(const char *uplo, mpackint n, mpackint nb, mpreal * A, mpackint lda, mpreal * e, mpreal * tau, mpreal * w, mpackint ldw);
void Rsterf(mpackint n, mpreal * d, mpreal * e, mpackint *info);
void Rorgtr(const char *uplo, mpackint n, mpreal * a, mpackint lda, mpreal * tau,
    mpreal * work, mpackint lwork, mpackint *info);
void Rgetrf ( mpackint m, mpackint n, mpreal * A, mpackint lda, mpackint *ipiv, mpackint *info );
void Rgetrs ( const char *trans, mpackint n, mpackint nrhs, mpreal * A, mpackint lda, mpackint *ipiv, mpreal * B, mpackint ldb, mpackint *info );
void Rlaswp ( mpackint n, mpreal * A, mpackint lda, mpackint k1, mpackint k2, mpackint *ipiv, mpackint incx );
void Rgetf2 ( mpackint m, mpackint n, mpreal * A, mpackint lda, mpackint *ipiv, mpackint *info );
#endif
