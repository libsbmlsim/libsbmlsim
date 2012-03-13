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
#define MAX_MTHD_NAME_LENGTH 256

#define EPSIRON 1.0e-8
#define DOUBLE_EQ(x, v) (((v - EPSIRON) < x) && (x < (v + EPSIRON)))

/* methods for numerical integration
 * n1: n'th order explicit method
 * n0: n'th order implicit method
 */
/* Explicit methods */
#define MTHD_EULER 1
#define MTHD_ADAMS_BASHFORTH_1 1
#define MTHD_ADAMS_BASHFORTH_2 11
#define MTHD_ADAMS_BASHFORTH_3 21
#define MTHD_ADAMS_BASHFORTH_4 31
#define MTHD_RUNGE_KUTTA 41
/* Name for explicit methods */
#define MTHD_NAME_EULER "Euler"
#define MTHD_NAME_ADAMS_BASHFORTH_1 "1st order Adams-Bashforth"
#define MTHD_NAME_ADAMS_BASHFORTH_2 "2nd order Adams-Bashforth"
#define MTHD_NAME_ADAMS_BASHFORTH_3 "3rd order Adams-Bashforth"
#define MTHD_NAME_ADAMS_BASHFORTH_4 "4th order Adams-Bashforth"
#define MTHD_NAME_RUNGE_KUTTA "4th order Runge-Kutta"

/* Implicit methods */
#define MTHD_BACKWARD_EULER 0
#define MTHD_CRANK_NICOLSON 10
#define MTHD_ADAMS_MOULTON_2 10
#define MTHD_ADAMS_MOULTON_3 20
#define MTHD_ADAMS_MOULTON_4 30
#define MTHD_BACKWARD_DIFFERENTIATION_2 40
#define MTHD_BACKWARD_DIFFERENTIATION_3 50
#define MTHD_BACKWARD_DIFFERENTIATION_4 60
/* Name for implicit methods */
#define MTHD_NAME_BACKWARD_EULER "Backward-Euler"
#define MTHD_NAME_CRANK_NICOLSON "Crank-Nicolson"
#define MTHD_NAME_ADAMS_MOULTON_2 "2nd order Adams-Moulton"
#define MTHD_NAME_ADAMS_MOULTON_3 "3rd order Adams-Moulton"
#define MTHD_NAME_ADAMS_MOULTON_4 "4th order Adams-Moulton"
#define MTHD_NAME_BACKWARD_DIFFERENTIATION_2 "2nd order Backward Differentiation"
#define MTHD_NAME_BACKWARD_DIFFERENTIATION_3 "3rd order Backward Differentiation"
#define MTHD_NAME_BACKWARD_DIFFERENTIATION_4 "4th order Backward Differentiation"

#endif  /* LibSBMLSim_Common_h */
