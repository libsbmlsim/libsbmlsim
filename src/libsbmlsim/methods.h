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
#define MTHD_RUNGE_KUTTA_FEHLBERG_5 51
#define MTHD_CASH_KARP 61


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
#define MTHD_NAME_RUNGE_KUTTA_FEHLBERG_5 "5th order Runge-Kutta-Fehlberg"
#define MTHD_NAME_CASH_KARP "5th order Cash-Karp"

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
