/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: mpack_config.h,v 1.7 2009/09/24 07:25:57 nakatamaho Exp $ 
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

/* work in progress */
/* put some definitons on mpack */

/* should depend on C compiler and environment 
   our intention is that use 64bit int when USE64BITINT is set.
   This should be the default on 64bit environment.
*/

#ifndef _MPACK_CONFIG_H_
#define _MPACK_CONFIG_H_

#include <stdlib.h>
#include <inttypes.h>

// #define USE64BITINT

#ifdef USE64BITINT
typedef int64_t mpackint;
#else
typedef int32_t mpackint;
#endif

#ifdef USE64BITINT
inline mpackint mpackabs(mpackint i)
{
  return labs(i);
}
#else
inline mpackint mpackabs(mpackint i)
{
  return abs(i);
}
#endif

typedef mpackint mpacklogical;

#ifdef __cplusplus
typedef mpacklogical(*ML_fp) (...);
#else
typedef mpacklogical(*ML_fp);
#endif

#endif
