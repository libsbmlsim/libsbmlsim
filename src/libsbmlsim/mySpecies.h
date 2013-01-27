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
#ifndef LibSBMLSim_MySpecies_h
#define LibSBMLSim_MySpecies_h

#include "typedefs.h"
#include "boolean.h"
#include "myRule.h"
#include "myCompartment.h"
#include <sbml/SBMLTypes.h>

struct _mySpecies {
  Species_t *origin;
  double value;
  double temp_value;
  boolean is_amount;
  boolean is_concentration;
  int has_only_substance_units;
  myCompartment *locating_compartment;
  double k[6]; /* for runge kutta */
  double **delay_val;
  unsigned int delay_val_width;
  unsigned int delay_val_length;
  myRule *depending_rule;
  double prev_val[3]; /* previous values for multistep solution */
  double prev_k[3]; /* previous values for multistep solution */
};

mySpecies *mySpecies_create();
void mySpecies_initWithModel(mySpecies *species, Model_t *model, int index);
void mySpecies_initDelayVal(mySpecies *species, unsigned int length, unsigned int width);
void mySpecies_free(mySpecies *species);

void mySpecies_reallocDelayVal(mySpecies *species, unsigned int length, unsigned int width);
void mySpecies_setLocatingCompartment(mySpecies *species, myCompartment *compartment);
void mySpecies_setDependingRule(mySpecies *species, myRule *rule);

#endif /* LibSBMLSim_MySpecies_h */
