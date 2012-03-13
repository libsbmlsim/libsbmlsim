#ifndef LibSBMLSim_Common_h
#define LibSBMLSim_Common_h

/* debug print */
#ifdef DEBUG_PRINT
#define dbg_printf printf
#else
#define dbg_printf 1 ? (void) 0 : printf
#endif

/* progress print */
#ifdef PROGRESS_PRINT
#define prg_printf printf
#else
#define prg_printf 1 ? (void) 0 : printf
#endif

/* defined variables */
#define MAX_MATH_LENGTH 256
#define MAX_DELAY_REACTION_NUM 256
#define MAX_ARG_NUM 256
#define MAX_ALLOCATED_MEMORY 1024
#define MAX_COPIED_AST 1024
#define MAX_ALGEBRAIC_VARIABLES 1024
#define MAX_ALGEBRAIC_CONSTANTS 1024
#define MAX_IDENTICAL_EVENTS 256
#define MAX_EVENTASSIGNMENTS 256
#define MAX_TIME_VARIANT_ASSIGNMENT 1024
#define MAX_INCLUDING_SPECIES 256

#define EPSIRON 1.0e-8
#define DOUBLE_EQ(x, v) (((v - EPSIRON) < x) && (x < (v + EPSIRON)))

#endif  /* LibSBMLSim_Common_h */
