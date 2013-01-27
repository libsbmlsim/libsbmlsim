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
#include "libsbmlsim/myParameter.h"
#include <stdlib.h>
#include <sbml/SBMLTypes.h>

myParameter *myParameter_create() {
  myParameter *parameter = (myParameter *)malloc(sizeof(myParameter));
  parameter->origin = NULL;
  parameter->value = 0;
  parameter->temp_value = 0;
  parameter->k[0] = 0;
  parameter->k[1] = 0;
  parameter->k[2] = 0;
  parameter->k[3] = 0;
  parameter->k[4] = 0;
  parameter->k[5] = 0;
  parameter->delay_val = NULL;
  parameter->delay_val_width = 0;
  parameter->delay_val_length = 0;
  parameter->depending_rule = NULL;
  parameter->prev_val[0] = 0;
  parameter->prev_val[1] = 0;
  parameter->prev_val[2] = 0;
  parameter->prev_k[0] = 0;
  parameter->prev_k[1] = 0;
  parameter->prev_k[2] = 0;
  return parameter;
}

void myParameter_initWithModel(myParameter *parameter, Model_t *model, int index) {
  Parameter_t *origin = (Parameter_t *)ListOf_get(Model_getListOfParameters(model), index);

  parameter->origin = origin;
  if (Parameter_isSetValue(origin)) {
    parameter->value = Parameter_getValue(origin);
  } else {
    parameter->value = 0;
  }
  parameter->temp_value = parameter->value;
  parameter->prev_val[0] = parameter->value;
  parameter->prev_val[1] = parameter->value;
  parameter->prev_val[2] = parameter->value;
}

void myParameter_initDelayVal(myParameter *parameter, unsigned int length, unsigned int width) {
  unsigned int i;
  parameter->delay_val_width = width;
  parameter->delay_val_length = length;
  parameter->delay_val = (double **)malloc(sizeof(double *) * length);
  for (i = 0; i < length; i++) {
    parameter->delay_val[i] = (double *)malloc(sizeof(double) * width);
  }
}

void myParameter_free(myParameter *parameter) {
  unsigned int i;

  if (parameter == NULL) {
    return;
  }

  if (parameter->delay_val != NULL) {
    for (i = 0; i < parameter->delay_val_length; i++) {
      free(parameter->delay_val[i]);
    }
    free(parameter->delay_val);
  }
  free(parameter);
}

void myParameter_reallocDelayVal(myParameter *parameter, unsigned int length, unsigned int width) {
  unsigned int i;
  unsigned int old_width = parameter->delay_val_width;
  unsigned int old_length = parameter->delay_val_length;
  double **delay_val = (double **)realloc(parameter->delay_val, sizeof(double *) * length);

  parameter->delay_val_width = width;
  parameter->delay_val_length = length;
  for (i = 0; i < length; i++) {
    if (i < old_length && width != old_width) {
      delay_val[i] = (double *)realloc(parameter->delay_val[i], sizeof(double) * width);
    } else if (i >= old_length) {
      delay_val[i] = (double *)malloc(sizeof(double) * width);
    }
  }
  parameter->delay_val = delay_val;
}

Parameter_t *myParameter_getOrigin(myParameter *parameter) {
  return parameter->origin;
}

void myParameter_setDependingRule(myParameter *parameter, myRule *rule) {
  parameter->depending_rule = rule;
}

