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
  mySpeciesReference *ret = (mySpeciesReference *)malloc(sizeof(mySpeciesReference));
  return ret;
}

void mySpeciesReference_initWithOrigin(mySpeciesReference *ref, SpeciesReference_t *origin) {
  ref->origin = origin;
  ref->mySp = NULL;
  ref->delay_val = NULL;
  ref->depending_rule = NULL;
  ref->eq = NULL;
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
  if (ref->eq != NULL) {
    equation_free(ref->eq);
  }
  free(ref);
}

void mySpeciesReference_setSpecies(mySpeciesReference *ref, mySpecies *species) {
  ref->mySp = species;
}

