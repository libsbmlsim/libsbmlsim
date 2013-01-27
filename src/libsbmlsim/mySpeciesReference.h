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
#ifndef LibSBMLSim_MySpeciesReference_h
#define LibSBMLSim_MySpeciesReference_h

#include "typedefs.h"
#include "equation.h"
#include "mySpecies.h"
#include "myRule.h"
#include "myReaction.h"
#include <sbml/SBMLTypes.h>

struct _mySpeciesReference {
  mySpecies *mySp;
  SpeciesReference_t *origin;
  equation *eq; /* for l2v4 */
  double value;
  double temp_value;
  double k[6]; /* for runge kutta */
  double **delay_val;
  unsigned int delay_val_width;
  unsigned int delay_val_length;
  myRule *depending_rule;
  double prev_val[3];
  double prev_k[3];
};

mySpeciesReference *mySpeciesReference_create();
void mySpeciesReference_initWithOrigin(mySpeciesReference *ref, SpeciesReference_t *origin);
void mySpeciesReference_initAsProduct(mySpeciesReference *ref, myReaction *reaction, int index);
void mySpeciesReference_initAsReactant(mySpeciesReference *ref, myReaction *reaction, int index);
void mySpeciesReference_free(mySpeciesReference *ref);

void mySpeciesReference_setSpecies(mySpeciesReference *ref, mySpecies *species);
void mySpeciesReference_setDependingRule(mySpeciesReference *ref, myRule *rule);

#endif /* LibSBMLSim_MySpeciesReference_h */
