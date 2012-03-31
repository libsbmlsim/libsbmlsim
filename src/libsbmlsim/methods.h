#ifndef LibSBMLSim_Methods_h
#define LibSBMLSim_Methods_h

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

/* Implicit methods */
#define MTHD_BACKWARD_EULER 0
#define MTHD_CRANK_NICOLSON 10
#define MTHD_ADAMS_MOULTON_2 10
#define MTHD_ADAMS_MOULTON_3 20
#define MTHD_ADAMS_MOULTON_4 30
#define MTHD_BACKWARD_DIFFERENCE_2 40
#define MTHD_BACKWARD_DIFFERENCE_3 50
#define MTHD_BACKWARD_DIFFERENCE_4 60

/* Name for explicit methods */
#define MTHD_NAME_EULER "Euler"
#define MTHD_NAME_ADAMS_BASHFORTH_1 "1st order Adams-Bashforth"
#define MTHD_NAME_ADAMS_BASHFORTH_2 "2nd order Adams-Bashforth"
#define MTHD_NAME_ADAMS_BASHFORTH_3 "3rd order Adams-Bashforth"
#define MTHD_NAME_ADAMS_BASHFORTH_4 "4th order Adams-Bashforth"
#define MTHD_NAME_RUNGE_KUTTA "4th order Runge-Kutta"

/* Name for implicit methods */
#define MTHD_NAME_BACKWARD_EULER "Backward-Euler"
#define MTHD_NAME_CRANK_NICOLSON "Crank-Nicolson"
#define MTHD_NAME_ADAMS_MOULTON_2 "2nd order Adams-Moulton"
#define MTHD_NAME_ADAMS_MOULTON_3 "3rd order Adams-Moulton"
#define MTHD_NAME_ADAMS_MOULTON_4 "4th order Adams-Moulton"
#define MTHD_NAME_BACKWARD_DIFFERENCE_2 "2nd order Backward Difference"
#define MTHD_NAME_BACKWARD_DIFFERENCE_3 "3rd order Backward Difference"
#define MTHD_NAME_BACKWARD_DIFFERENCE_4 "4th order Backward Difference"

#endif  /* LibSBMLSim_Methods_h */
