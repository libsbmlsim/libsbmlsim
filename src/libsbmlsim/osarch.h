#ifndef LibSBMLSim_OsArch_h
#define LibSBMLSim_OsArch_h

/* handle int64_t */
#if defined(_MSC_VER) && (_MSC_VER <= 1500)
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;
#else
#include <stdint.h>
#endif

#endif  /* LibSBMLSim_OsArch_h */
