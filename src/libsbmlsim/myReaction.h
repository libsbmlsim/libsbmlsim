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
#ifndef LibSBMLSim_MyReaction_h
#define LibSBMLSim_MyReaction_h

#include "typedefs.h"
#include "boolean.h"
#include "mySpeciesReference.h"
#include <sbml/SBMLTypes.h>

struct _myReaction {
  Reaction_t *origin;
  equation *eq;
  mySpeciesReference **products;
  unsigned int num_of_products;
  mySpeciesReference **reactants;
  unsigned int num_of_reactants;
  boolean is_fast;
  boolean is_reversible;
  equation *products_equili_numerator;
  equation *reactants_equili_numerator;
};

myReaction *myReaction_create();
void myReaction_initWithModel(myReaction *reaction, Model_t *model, int index);
void myReaction_free(myReaction *reaction);

void myReaction_initProducts(myReaction *reaction, mySpecies **species, unsigned int num_of_species);
void myReaction_initReactants(myReaction *reaction, mySpecies **species, unsigned int num_of_species);

Reaction_t *myReaction_getOrigin(myReaction *reaction);
void myReaction_addProduct(myReaction *reaction, mySpeciesReference *product);
void myReaction_addReactant(myReaction *reaction, mySpeciesReference *reactant);

#endif /* LibSBMLSim_MyReaction_h */
