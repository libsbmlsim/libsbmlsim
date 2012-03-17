#ifndef LibSBMLSim_OsArch_h
#define LibSBMLSim_OsArch_h

#include <sys/types.h>

/* handle int64_t */
#if defined(_MSC_VER) && (_MSC_VER <= 1500)
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;
#else
#include <stdint.h>
#endif

/* Endian for Windows */
/* (IA-64 is bi-endian, but Windows use it in little endian) */
#ifdef _WIN32
#define _LITTLE_ENDIAN    1234
#define BYTE_ORDER _LITTLE_ENDIAN
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

#endif  /* LibSBMLSim_OsArch_h */
