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
#include "libsbmlsim/myInitialAssignment.h"
#include <stdlib.h>
#include <string.h>
#include <sbml/SBMLTypes.h>

myInitialAssignment *myInitialAssignment_create() {
  myInitialAssignment *assign = (myInitialAssignment *)malloc(sizeof(myInitialAssignment));
  assign->origin = NULL;
  assign->eq = NULL;
  assign->target_species = NULL;
  assign->target_parameter = NULL;
  assign->target_compartment = NULL;
  assign->target_species_reference = NULL;
  return assign;
}

void myInitialAssignment_initWithModel(myInitialAssignment *assign, Model_t *model, int index) {
  InitialAssignment_t *origin = (InitialAssignment_t*)ListOf_get(Model_getListOfInitialAssignments(model), index);

  assign->origin = origin;
  assign->eq = equation_create();
}

void myInitialAssignment_free(myInitialAssignment *assign) {
  if (assign == NULL) {
    return;
  }

  if (assign->eq != NULL) {
    equation_free(assign->eq);
  }
  free(assign);
}

void myInitialAssignment_initTarget(
    myInitialAssignment *assign,
    mySpecies **species, unsigned int num_of_species,
    myParameter **parameters, unsigned int num_of_parameters,
    myCompartment **compartments, unsigned int num_of_compartments) {
  unsigned int i;
  InitialAssignment_t *origin;
  const char *origin_symbol;

  origin = myInitialAssignment_getOrigin(assign);
  origin_symbol = InitialAssignment_getSymbol(origin);

  for (i = 0; i < num_of_species; i++) {
    if (strcmp(origin_symbol, Species_getId(species[i]->origin)) == 0) {
      myInitialAssignment_setTargetSpecies(assign, species[i]);
      return;
    }
  }

  for (i = 0; i < num_of_parameters; i++) {
    if (strcmp(origin_symbol, Parameter_getId(parameters[i]->origin)) == 0) {
      myInitialAssignment_setTargetParameter(assign, parameters[i]);
      return;
    }
  }

  for (i = 0; i < num_of_compartments; i++) {
    if (strcmp(origin_symbol, Compartment_getId(compartments[i]->origin)) == 0) {
      myInitialAssignment_setTargetCompartment(assign, compartments[i]);
      return;
    }
  }
}

void myInitialAssignment_initTargetSpeciesReference(
    myInitialAssignment *assign,
    myReaction **reactions, unsigned int num_of_reactions) {
  unsigned int i, j;
  InitialAssignment_t *origin;
  myReaction *re;
  mySpeciesReference *product, *reactant;
  const char *origin_symbol;

  origin = myInitialAssignment_getOrigin(assign);
  origin_symbol = InitialAssignment_getSymbol(origin);

  for (i = 0; i < num_of_reactions; i++) {
    re = reactions[i];

    for (j = 0; j < re->num_of_products; j++) {
      product = re->products[j];
      if (SpeciesReference_isSetId(product->origin)
          && strcmp(origin_symbol, SpeciesReference_getId(product->origin)) == 0) {
        myInitialAssignment_setTargetSpeciesReference(assign, product);
        return;
      }
    }

    for (j = 0; j < re->num_of_reactants; j++) {
      reactant = re->reactants[j];
      if (SpeciesReference_isSetId(reactant->origin)
          && strcmp(origin_symbol, SpeciesReference_getId(reactant->origin)) == 0) {
        myInitialAssignment_setTargetSpeciesReference(assign, reactant);
        return;
      }
    }
  }
}

InitialAssignment_t *myInitialAssignment_getOrigin(myInitialAssignment *assign) {
  return assign->origin;
}

void myInitialAssignment_setTargetSpecies(myInitialAssignment *assign, mySpecies *species) {
  assign->target_species = species;
}

void myInitialAssignment_setTargetParameter(myInitialAssignment *assign, myParameter *parameter) {
  assign->target_parameter = parameter;
}

void myInitialAssignment_setTargetCompartment(myInitialAssignment *assign, myCompartment *compartment) {
  assign->target_compartment = compartment;
}

void myInitialAssignment_setTargetSpeciesReference(myInitialAssignment *assign, mySpeciesReference *ref) {
  assign->target_species_reference = ref;
}

void myInitialAssignment_setEquation(myInitialAssignment *assign, equation *eq) {
  assign->eq = eq;
}

