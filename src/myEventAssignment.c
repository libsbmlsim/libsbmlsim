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
#include "libsbmlsim/myEventAssignment.h"
#include <stdlib.h>
#include <string.h>
#include <sbml/SBMLTypes.h>

myEventAssignment *myEventAssignment_create() {
  myEventAssignment *assign = (myEventAssignment *)malloc(sizeof(myEventAssignment));
  assign->origin = NULL;
  assign->eq = NULL;
  assign->target_species = NULL;
  assign->target_parameter = NULL;
  assign->target_compartment = NULL;
  assign->target_species_reference = NULL;
  return assign;
}

void myEventAssignment_initWithEvent(myEventAssignment *assign, Event_t *event, int index) {
  EventAssignment_t *origin = (EventAssignment_t *)ListOf_get(Event_getListOfEventAssignments(event), index);

  assign->origin = origin;
  assign->eq = equation_create();
}

void myEventAssignment_free(myEventAssignment *assign) {
  if (assign == NULL) {
    return;
  }
  
  if (assign->eq != NULL) {
    equation_free(assign->eq);
  }
  free(assign);
}

void myEventAssignment_initTarget(
    myEventAssignment *assign,
    mySpecies **species, unsigned int num_of_species,
    myParameter **parameters, unsigned int num_of_parameters,
    myCompartment **compartments, unsigned int num_of_compartments,
    myReaction **reactions, unsigned int num_of_reactions) {
  unsigned int i, j;
  EventAssignment_t *origin;
  const char *origin_var;
  mySpecies *sp;
  myParameter *param;
  myCompartment *comp;
  myReaction *re;
  mySpeciesReference *product, *reactant;

  origin = myEventAssignment_getOrigin(assign);
  origin_var = EventAssignment_getVariable(origin);

  for (i = 0; i < num_of_species; i++) {
    sp = species[i];
    if (strcmp(origin_var, Species_getId(sp->origin)) == 0) {
      myEventAssignment_setTargetSpecies(assign, sp);
      return;
    }
  }

  for (i = 0; i < num_of_parameters; i++) {
    param = parameters[i];
    if (strcmp(origin_var, Parameter_getId(param->origin)) == 0) {
      myEventAssignment_setTargetParameter(assign, param);
      return;
    }
  }

  for (i = 0; i < num_of_compartments; i++) {
    comp = compartments[i];
    if (strcmp(origin_var, Compartment_getId(comp->origin)) == 0) {
      myEventAssignment_setTargetCompartment(assign, comp);
      return;
    }
  }

  for (i = 0; i < num_of_reactions; i++) {
    re = reactions[i];

    for (j = 0; j < re->num_of_products; j++) {
      product = re->products[j];
      if (SpeciesReference_isSetId(product->origin)
          && strcmp(origin_var, SpeciesReference_getId(product->origin)) == 0) {
        myEventAssignment_setTargetSpeciesReference(assign, product);
        return;
      }
    }

    for (j = 0; j < re->num_of_reactants; j++) {
      reactant = re->reactants[j];
      if (SpeciesReference_isSetId(reactant->origin)
          && strcmp(origin_var, SpeciesReference_getId(reactant->origin)) == 0) {
        myEventAssignment_setTargetSpeciesReference(assign, reactant);
        return;
      }
    }
  }
}

EventAssignment_t *myEventAssignment_getOrigin(myEventAssignment *assign) {
  return assign->origin;
}

void myEventAssignment_setTargetSpecies(myEventAssignment *assign, mySpecies *species) {
  assign->target_species = species;
}

void myEventAssignment_setTargetParameter(myEventAssignment *assign, myParameter *parameter) {
  assign->target_parameter = parameter;
}

void myEventAssignment_setTargetCompartment(myEventAssignment *assign, myCompartment *compartment) {
  assign->target_compartment = compartment;
}

void myEventAssignment_setTargetSpeciesReference(myEventAssignment *assign, mySpeciesReference *ref) {
  assign->target_species_reference = ref;
}

