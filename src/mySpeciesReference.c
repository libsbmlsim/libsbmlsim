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
#include "libsbmlsim/mySpeciesReference.h"
#include <stdlib.h>
#include <sbml/SBMLTypes.h>

mySpeciesReference *mySpeciesReference_create() {
  mySpeciesReference *ref = (mySpeciesReference *)malloc(sizeof(mySpeciesReference));
  ref->mySp = NULL;
  ref->origin = NULL;
  ref->eq = NULL;
  ref->value = 0;
  ref->temp_value = 0;
  ref->k[0] = 0;
  ref->k[1] = 0;
  ref->k[2] = 0;
  ref->k[3] = 0;
  ref->k[4] = 0;
  ref->k[5] = 0;
  ref->delay_val = NULL;
  ref->delay_val_width = 0;
  ref->delay_val_length = 0;
  ref->depending_rule = NULL;
  ref->prev_val[0] = 0;
  ref->prev_val[1] = 0;
  ref->prev_val[2] = 0;
  ref->prev_k[0] = 0;
  ref->prev_k[1] = 0;
  ref->prev_k[2] = 0;
  return ref;
}

void mySpeciesReference_initWithOrigin(mySpeciesReference *ref, SpeciesReference_t *origin) {
  ref->origin = origin;
  ref->eq = equation_create();
}

void mySpeciesReference_initAsReactant(mySpeciesReference *ref, myReaction *reaction, int index) {
  SpeciesReference_t *origin = (SpeciesReference_t *)Reaction_getReactant(reaction->origin, index);
  mySpeciesReference_initWithOrigin(ref, origin);
}

void mySpeciesReference_initAsProduct(mySpeciesReference *ref, myReaction *reaction, int index) {
  SpeciesReference_t *origin = (SpeciesReference_t *)Reaction_getProduct(reaction->origin, index);
  mySpeciesReference_initWithOrigin(ref, origin);
}

void mySpeciesReference_free(mySpeciesReference *ref) {
  if (ref == NULL) {
    return;
  }

  if (ref->eq != NULL) {
    equation_free(ref->eq);
  }
  free(ref);
}

void mySpeciesReference_setSpecies(mySpeciesReference *ref, mySpecies *species) {
  ref->mySp = species;
}

void mySpeciesReference_setDependingRule(mySpeciesReference *ref, myRule *rule) {
  ref->depending_rule = rule;
}

