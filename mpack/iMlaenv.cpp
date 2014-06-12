/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: iMlaenv.cpp,v 1.5 2009/09/12 07:59:10 nakatamaho Exp $ 
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
#include <string.h>
#include <ctype.h>

#define MLANAMESIZE 6

//ISPEC = 1:  block size
//In these examples, separate code is provided for setting NB for
//real and complex.  We assume that NB will take the same value in
//single or double precision.

mpackint
iMlaenv1(const char *Mlaname, const char *opts, mpackint n1, mpackint n2, mpackint n3, mpackint n4)
{
    mpackint nb = 1;
#if !defined (IMLAENV_DEBUG)
    if (strcmp(&Mlaname[1],"orgqr") == 0) { nb = 32; return nb; }
    if (strcmp(&Mlaname[1],"orgql") == 0) { nb = 32; return nb; }
    if (strcmp(&Mlaname[1],"potrf") == 0) { nb = 64; return nb; }
    if (strcmp(&Mlaname[1],"trtri") == 0) { nb = 64; return nb; }
    if (strcmp(&Mlaname[1],"dsytrd") == 0) { nb = 32;return nb;  }
    if (strcmp(&Mlaname[1],"getrf") == 0)  { nb = 64;return nb;  }
    if (strcmp(&Mlaname[1],"getri") == 0)  { nb = 64;return nb;  }
#else
    if (strcmp(&Mlaname[1],"potrf") == 0)  { nb = 8;return nb; }
    if (strcmp(&Mlaname[1],"orgqr") == 0)  { nb = 8;return nb; }
    if (strcmp(&Mlaname[1],"orgql") == 0)  { nb = 8;return nb; }
    if (strcmp(&Mlaname[1],"trtri") == 0)  { nb = 8;return nb; }
    if (strcmp(&Mlaname[0],"dsytrd") == 0) { nb = 8;return nb; }
    if (strcmp(&Mlaname[1],"getrf") == 0)  { nb = 8;return nb; }
    if (strcmp(&Mlaname[1],"getri") == 0)  { nb = 8;return nb; }
#endif
    return nb;
}

//*     ISPEC = 2: minimum block size
mpackint
iMlaenv2(const char *Mlaname, const char *opts, mpackint n1, mpackint n2, mpackint n3, mpackint n4)
{
    mpackint nbmin = 1;
    if (strcmp(&Mlaname[1], "orgqr") == 0)  { nbmin = 2; return nbmin; }
    if (strcmp(&Mlaname[1], "orgql") == 0)  { nbmin = 2; return nbmin; }
    if (strcmp(&Mlaname[1], "trtri") == 0)  { nbmin = 2; return nbmin; }
    if (strcmp(&Mlaname[0], "dsytrd") == 0) { nbmin = 2; return nbmin; }
    if (strcmp(&Mlaname[0], "getri") == 0)  { nbmin = 2; return nbmin; }

    return nbmin;
}

//     ISPEC = 3:  crossover point
mpackint
iMlaenv3(const char *Mlaname, const char *opts, mpackint n1, mpackint n2, mpackint n3, mpackint n4)
{
    mpackint nx = 1;
#if !defined (IMLAENV_DEBUG)
    if (strcmp(&Mlaname[1],"orgqr")==0) { nx = 128; return nx; }
    if (strcmp(&Mlaname[1],"orgql")==0) { nx = 128; return nx; }
    if (strcmp(&Mlaname[0],"dsytrd")==0){ nx = 32; return nx; }
#else
    if (strcmp(&Mlaname[1],"orgqr") == 0) { nx = 6; return nx; }
    if (strcmp(&Mlaname[1],"orgql") == 0) { nx = 6; return nx; }
    if (strcmp(&Mlaname[0], "dsytrd")== 0){ nx = 6; return nx; }
#endif
    return nx;
}

mpackint
iMlaenv4(const char *Mlaname, const char *opts, mpackint n1, mpackint n2, mpackint n3, mpackint n4)
{
    return 1;
}

mpackint
iMlaenv5(const char *Mlaname, const char *opts, mpackint n1, mpackint n2, mpackint n3, mpackint n4)
{
    return 1;
}

mpackint
iMlaenv6(const char *Mlaname, const char *opts, mpackint n1, mpackint n2, mpackint n3, mpackint n4)
{
    return 1;
}

mpackint
iMlaenv7(const char *Mlaname, const char *opts, mpackint n1, mpackint n2, mpackint n3, mpackint n4)
{
    return 1;
}

mpackint
iMlaenv8(const char *Mlaname, const char *opts, mpackint n1, mpackint n2, mpackint n3, mpackint n4)
{
    return 1;
}

mpackint
iMlaenv9(const char *Mlaname, const char *opts, mpackint n1, mpackint n2, mpackint n3, mpackint n4)
{
    return 1;
}

mpackint
iMlaenv10(const char *Mlaname, const char *opts, mpackint n1, mpackint n2, mpackint n3,
    mpackint n4)
{
    return 1;
}

mpackint
iMlaenv11(const char *Mlaname, const char *opts, mpackint n1, mpackint n2, mpackint n3,
    mpackint n4)
{
    return 1;
}

mpackint
iMlaenv12(const char *Mlaname, const char *opts, mpackint n1, mpackint n2, mpackint n3,
    mpackint n4)
{
    return 1;
}

mpackint
iMlaenv13(const char *Mlaname, const char *opts, mpackint n1, mpackint n2, mpackint n3,
    mpackint n4)
{
    return 1;
}

mpackint
iMlaenv14(const char *Mlaname, const char *opts, mpackint n1, mpackint n2, mpackint n3,
    mpackint n4)
{
    return 1;
}

mpackint
iMlaenv15(const char *Mlaname, const char *opts, mpackint n1, mpackint n2, mpackint n3,
    mpackint n4)
{
    return 1;
}

mpackint
iMlaenv16(const char *Mlaname, const char *opts, mpackint n1, mpackint n2, mpackint n3,
    mpackint n4)
{
    return 1;
}

mpackint
iMlaenv_gmp(mpackint ispec, const char *name, const char *opts, mpackint n1, mpackint n2, mpackint n3,
    mpackint n4)
{
    mpackint iret, i, up;

    iret = -1;

    char Mlaname[MLANAMESIZE + 1];
//buggy
    strncpy(Mlaname, name, MLANAMESIZE);
    for (i = 0; i < MLANAMESIZE; i++) {
	up = tolower(Mlaname[i]);
	Mlaname[i] = up;
    }

    if (!Mlsame_gmp(Mlaname, "r") && !Mlsame_gmp(Mlaname, "c"))
	return iret;

    switch (ispec) {
    case 1:
	iret = iMlaenv1(Mlaname, opts, n1, n2, n3, n4);
	break;
    case 2:
	iret = iMlaenv2(Mlaname, opts, n1, n2, n3, n4);
	break;
    case 3:
	iret = iMlaenv3(Mlaname, opts, n1, n2, n3, n4);
	break;
    case 4:
	iret = iMlaenv4(Mlaname, opts, n1, n2, n3, n4);
	break;
    case 5:
	iret = iMlaenv5(Mlaname, opts, n1, n2, n3, n4);
	break;
    case 6:
	iret = iMlaenv6(Mlaname, opts, n1, n2, n3, n4);
	break;
    case 7:
	iret = iMlaenv7(Mlaname, opts, n1, n2, n3, n4);
	break;
    case 8:
	iret = iMlaenv8(Mlaname, opts, n1, n2, n3, n4);
	break;
    case 9:
	iret = iMlaenv9(Mlaname, opts, n1, n2, n3, n4);
	break;
    case 10:
	iret = iMlaenv10(Mlaname, opts, n1, n2, n3, n4);
	break;
    case 11:
	iret = iMlaenv11(Mlaname, opts, n1, n2, n3, n4);
	break;
    case 12:
	iret = iMlaenv12(Mlaname, opts, n1, n2, n3, n4);
	break;
    case 13:
	iret = iMlaenv13(Mlaname, opts, n1, n2, n3, n4);
	break;
    case 14:
	iret = iMlaenv14(Mlaname, opts, n1, n2, n3, n4);
	break;
    case 15:
	iret = iMlaenv15(Mlaname, opts, n1, n2, n3, n4);
	break;
    case 16:
	iret = iMlaenv16(Mlaname, opts, n1, n2, n3, n4);
	break;
    default:
	iret = -1;
    }
    return iret;
}
