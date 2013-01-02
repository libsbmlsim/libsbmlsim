/**
 * <!--------------------------------------------------------------------------
 * This file is part of libSBMLSim.  Please visit
 * http://fun.bio.keio.ac.jp/software/libsbmlsim/ for more
 * information about libSBMLSim and its latest version.
 *
 * Copyright (C) 2011-2013 by the Keio University, Yokohama, Japan
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution.
 * ---------------------------------------------------------------------- -->*/
#ifndef LibSBMLSim_OsArch_h
#define LibSBMLSim_OsArch_h

#include <sys/types.h>

/* handle int64_t, u_int32_t etc. */
#if defined(_MSC_VER) && (_MSC_VER <= 1500)
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;
#else
#include <stdint.h>
#endif

/* Endian for Windows */
/* (IA-64 is bi-endian, but Windows use it in little endian) */
#ifdef _WIN32
#define LITTLE_ENDIAN    1234
#define BIG_ENDIAN    4321
#define BYTE_ORDER LITTLE_ENDIAN
#elif defined(__linux__)
#include <endian.h>
#else    /* __linux__ */
#include <machine/endian.h>
#endif

/* asinh(), acosh(), atanh(), log1p is not C90 standard */
#if defined(_MSC_VER) || defined(__STRICT_ANSI__)
#define __USE_BSD
#include "math_private.h"
#endif

/* DLL export */
#ifdef _WIN32
#define SBMLSIM_EXPORT __declspec(dllexport)
#else
#define SBMLSIM_EXPORT
#endif

#endif  /* LibSBMLSim_OsArch_h */
