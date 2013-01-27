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
#ifndef LibSBMLSim_MyCompartment_h
#define LibSBMLSim_MyCompartment_h

#include "typedefs.h"
#include "common.h"
#include "mySpecies.h"
#include "myRule.h"
#include <sbml/SBMLTypes.h>

struct _myCompartment {
  Compartment_t* origin;
  double value; /* compartment "size" value */
  double temp_value;
  double k[6]; /* for runge kutta */
  double **delay_val;
  unsigned int delay_val_width;
  unsigned int delay_val_length;
  myRule *depending_rule;
  double prev_val[3]; /* previous values for multistep solution */
  double prev_k[3]; /* previous values for multistep solution */
  mySpecies *including_species[MAX_INCLUDING_SPECIES];
  unsigned int num_of_including_species;
};

myCompartment *myCompartment_create();
void myCompartment_initWithModel(myCompartment *compartment, Model_t *model, int index);
void myCompartment_initDelayVal(myCompartment *compartment, unsigned int length, unsigned int width);
void myCompartment_free(myCompartment *compartment);

void myCompartment_reallocDelayVal(myCompartment *compartment, unsigned int length, unsigned int width);
Compartment_t *myCompartment_getOrigin(myCompartment *compartment);
void myCompartment_setDependingRule(myCompartment *compartment, myRule *rule);
void myCompartment_addIncludingSpecies(myCompartment *compartment, mySpecies *species);

#endif /* LibSBMLSim_MyCompartment_h */
