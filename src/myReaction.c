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
#include "libsbmlsim/myReaction.h"
#include <stdlib.h>
#include <sbml/SBMLTypes.h>

myReaction *myReaction_create() {
  myReaction *ret = (myReaction *)malloc(sizeof(myReaction));
  return ret;
}

void myReaction_initWithModel(myReaction *reaction, Model_t *model, int index) {
  int num_products, num_reactants;
  Reaction_t *origin = (Reaction_t *)ListOf_get(Model_getListOfReactions(model), index);

  num_products = Reaction_getNumProducts(origin);
  num_reactants = Reaction_getNumReactants(origin);

  reaction->origin = origin;
  reaction->is_fast = Reaction_getFast(origin);
  reaction->is_reversible = Reaction_getReversible(origin);
  reaction->num_of_products = 0;
  reaction->num_of_reactants = 0;
  reaction->products = (mySpeciesReference **)malloc(sizeof(mySpeciesReference *) * num_products);
  reaction->reactants = (mySpeciesReference **)malloc(sizeof(mySpeciesReference *) * num_reactants);
  reaction->products_equili_numerator = NULL;
  reaction->reactants_equili_numerator = NULL;
  reaction->eq = NULL;
}

void myReaction_free(myReaction *reaction) {
  unsigned int i;

  for (i = 0; i < reaction->num_of_products; i++) {
    mySpeciesReference_free(reaction->products[i]);
  }
  free(reaction->products);

  for (i = 0; i < reaction->num_of_reactants; i++) {
    mySpeciesReference_free(reaction->reactants[i]);
  }
  free(reaction->reactants);

  if (reaction->products_equili_numerator != NULL) {
    equation_free(reaction->products_equili_numerator);
  }
  if (reaction->reactants_equili_numerator != NULL) {
    equation_free(reaction->reactants_equili_numerator);
  }
  if (reaction->eq != NULL) {
    equation_free(reaction->eq);
  }

  free(reaction);
}

void myReaction_addProduct(myReaction *reaction, mySpeciesReference *product) {
  reaction->products[reaction->num_of_products] = product;
  reaction->num_of_products++;
}

void myReaction_addReactant(myReaction *reaction, mySpeciesReference *reactant) {
  reaction->reactants[reaction->num_of_reactants] = reactant;
  reaction->num_of_reactants++;
}

