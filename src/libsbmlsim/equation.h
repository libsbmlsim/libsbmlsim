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
#ifndef LibSBMLSim_Equation_h
#define LibSBMLSim_Equation_h

#include "typedefs.h"
#include "common.h"
#include "boolean.h"

struct _equation {
  double *number[MAX_MATH_LENGTH];
  int op[MAX_MATH_LENGTH];
  double **delay_number[MAX_MATH_LENGTH];
  double **delay_comp_size[MAX_MATH_LENGTH];
  equation *explicit_delay_eq[MAX_MATH_LENGTH];
  unsigned int math_length;
  /* new code */
  boolean time_reverse_flag;
  double *reverse_time;
  /* new code end */
};

equation *equation_create();
void equation_free(equation *eq);

#endif /* LibSBMLSim_Equation_h */
