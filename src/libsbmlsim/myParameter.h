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
#ifndef LibSBMLSim_MyParameter_h
#define LibSBMLSim_MyParameter_h

#include "typedefs.h"
#include "myRule.h"
#include <sbml/SBMLTypes.h>

struct _myParameter {
  Parameter_t* origin;
  double value;
  double temp_value;
  double k[6]; /* for runge kutta */
  double **delay_val;
  unsigned int delay_val_width;
  unsigned int delay_val_length;
  myRule *depending_rule;
  double prev_val[3]; /* previous values for multistep solution */
  double prev_k[3]; /* previous values for multistep solution */
};

myParameter *myParameter_create();
void myParameter_initWithModel(myParameter *parameter, Model_t *model, int index);
void myParameter_initDelayVal(myParameter *parameter, unsigned int length, unsigned int width);
void myParameter_free(myParameter *parameter);

void myParameter_reallocDelayVal(myParameter *parameter, unsigned int length, unsigned int width);
Parameter_t *myParameter_getOrigin(myParameter *parameter);
void myParameter_setDependingRule(myParameter *parameter, myRule *rule);

#endif /* LibSBMLSim_MyParameter_h */
