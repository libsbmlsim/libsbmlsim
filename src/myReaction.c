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
#include <string.h>
#include <sbml/SBMLTypes.h>

myReaction *myReaction_create() {
  myReaction *reaction = (myReaction *)malloc(sizeof(myReaction));
  reaction->origin = NULL;
  reaction->eq = NULL;
  reaction->products = NULL;
  reaction->num_of_products = 0;
  reaction->reactants = NULL;
  reaction->num_of_reactants = 0;
  reaction->is_fast = false;
  reaction->is_reversible = false;
  reaction->products_equili_numerator = NULL;
  reaction->reactants_equili_numerator = NULL;
  return reaction;
}

void myReaction_initWithModel(myReaction *reaction, Model_t *model, int index) {
  int num_of_products, num_of_reactants;
  Reaction_t *origin = (Reaction_t *)ListOf_get(Model_getListOfReactions(model), index);

  num_of_products = Reaction_getNumProducts(origin);
  num_of_reactants = Reaction_getNumReactants(origin);

  reaction->origin = origin;
  reaction->eq = equation_create();
  reaction->products = (mySpeciesReference **)malloc(sizeof(mySpeciesReference *) * num_of_products);
  reaction->reactants = (mySpeciesReference **)malloc(sizeof(mySpeciesReference *) * num_of_reactants);
  reaction->is_fast = Reaction_getFast(origin);
  reaction->is_reversible = Reaction_getReversible(origin);
}

void myReaction_free(myReaction *reaction) {
  unsigned int i;

  if (reaction == NULL) {
    return;
  }

  if (reaction->eq != NULL) {
    equation_free(reaction->eq);
  }
  if (reaction->products != NULL) {
    for (i = 0; i < reaction->num_of_products; i++) {
      mySpeciesReference_free(reaction->products[i]);
    }
    free(reaction->products);
  }
  if (reaction->reactants != NULL) {
    for (i = 0; i < reaction->num_of_reactants; i++) {
      mySpeciesReference_free(reaction->reactants[i]);
    }
    free(reaction->reactants);
  }
  if (reaction->products_equili_numerator != NULL) {
    equation_free(reaction->products_equili_numerator);
  }
  if (reaction->reactants_equili_numerator != NULL) {
    equation_free(reaction->reactants_equili_numerator);
  }

  free(reaction);
}

void myReaction_initProducts(myReaction *reaction, mySpecies **species, unsigned int num_of_species) {
  unsigned int i, j;
  unsigned int num_of_products;
  Reaction_t *origin;
  const char *product_id;
  mySpeciesReference *ref;

  origin = myReaction_getOrigin(reaction);
  num_of_products = Reaction_getNumProducts(origin);

  for (i = 0; i < num_of_products; i++) {
    product_id = SpeciesReference_getSpecies(Reaction_getProduct(origin, i));
    for (j = 0; j < num_of_species; j++) {
      if (strcmp(product_id, Species_getId(species[j]->origin)) == 0) {
        ref = mySpeciesReference_create();
        mySpeciesReference_initAsProduct(ref, reaction, i);
        mySpeciesReference_setSpecies(ref, species[j]);
        myReaction_addProduct(reaction, ref);
      }
    }
  }
}

void myReaction_initReactants(myReaction *reaction, mySpecies **species, unsigned int num_of_species) {
  unsigned int i, j;
  unsigned int num_of_reactants;
  Reaction_t *origin;
  const char *reactant_id;
  mySpeciesReference *ref;

  origin = myReaction_getOrigin(reaction);
  num_of_reactants = Reaction_getNumReactants(origin);

  for (i = 0; i < num_of_reactants; i++) {
    reactant_id = SpeciesReference_getSpecies(Reaction_getReactant(origin, i));
    for (j = 0; j < num_of_species; j++) {
      if (strcmp(reactant_id, Species_getId(species[j]->origin)) == 0) {
        ref = mySpeciesReference_create();
        mySpeciesReference_initAsReactant(ref, reaction, i);
        mySpeciesReference_setSpecies(ref, species[j]);
        myReaction_addReactant(reaction, ref);
      }
    }
  }
}

Reaction_t *myReaction_getOrigin(myReaction *reaction) {
  return reaction->origin;
}

void myReaction_addProduct(myReaction *reaction, mySpeciesReference *product) {
  reaction->products[reaction->num_of_products] = product;
  reaction->num_of_products++;
}

void myReaction_addReactant(myReaction *reaction, mySpeciesReference *reactant) {
  reaction->reactants[reaction->num_of_reactants] = reactant;
  reaction->num_of_reactants++;
}

