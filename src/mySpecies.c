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
#include "libsbmlsim/mySpecies.h"
#include <stdlib.h>
#include <sbml/SBMLTypes.h>

mySpecies *mySpecies_create() {
  mySpecies *species = (mySpecies *)malloc(sizeof(mySpecies));
  species->origin = NULL;
  species->value = 0;
  species->temp_value = 0;
  species->is_amount = true;
  species->is_concentration = false;
  species->has_only_substance_units = false;
  species->locating_compartment = NULL;
  species->k[0] = 0;
  species->k[1] = 0;
  species->k[2] = 0;
  species->k[3] = 0;
  species->k[4] = 0;
  species->k[5] = 0;
  species->delay_val = NULL;
  species->delay_val_width = 0;
  species->delay_val_length = 0;
  species->depending_rule = NULL;
  species->prev_val[0] = 0;
  species->prev_val[1] = 0;
  species->prev_val[2] = 0;
  species->prev_k[0] = 0;
  species->prev_k[1] = 0;
  species->prev_k[2] = 0;
  return species;
}

void mySpecies_initWithModel(mySpecies *species, Model_t *model, int index) {
  Species_t *origin = (Species_t *)ListOf_get(Model_getListOfSpecies(model), index);

  species->origin = origin;

  if (Species_isSetInitialAmount(origin)) {
    species->value = Species_getInitialAmount(origin);
    species->is_amount = true;
    species->is_concentration = false;
  } else if (Species_isSetInitialConcentration(origin)) {
    species->value = Species_getInitialConcentration(origin);
    species->is_amount = false;
    species->is_concentration = true;
  } else if (Species_getHasOnlySubstanceUnits(origin) ||
      Compartment_getSpatialDimensions(Model_getCompartmentById(model, Species_getCompartment(origin))) == 0) {
    species->value = 0;
    species->is_amount = true;
    species->is_concentration = false;      
  } else {
    species->value = 0;
    species->is_amount = false;
    species->is_concentration = true;
  }

  species->has_only_substance_units = Species_getHasOnlySubstanceUnits(origin);
  species->temp_value = species->value;
  species->prev_val[0] = species->value;
  species->prev_val[1] = species->value;
  species->prev_val[2] = species->value;
}

void mySpecies_initDelayVal(mySpecies *species, unsigned int length, unsigned int width) {
  unsigned int i;
  species->delay_val_width = width;
  species->delay_val_length = length;
  species->delay_val = (double **)malloc(sizeof(double *) * length);
  for (i = 0; i < length; i++) {
    species->delay_val[i] = (double *)malloc(sizeof(double) * width);
  }
}

void mySpecies_free(mySpecies *species) {
  unsigned int i;

  if (species == NULL) {
    return;
  }

  if (species->delay_val != NULL) {
    for (i = 0; i < species->delay_val_length; i++) {
      free(species->delay_val[i]);
    }
    free(species->delay_val);
  }
  free(species);
}

void mySpecies_reallocDelayVal(mySpecies *species, unsigned int length, unsigned int width) {
  unsigned int i;
  unsigned int old_width = species->delay_val_width;
  unsigned int old_length = species->delay_val_length;
  double **delay_val = (double **)realloc(species->delay_val, sizeof(double *) * length);

  species->delay_val_width = width;
  species->delay_val_length = length;
  for (i = 0; i < length; i++) {
    if (i < old_length && width != old_width) {
      delay_val[i] = (double *)realloc(species->delay_val[i], sizeof(double) * width);
    } else if (i >= old_length) {
      delay_val[i] = (double *)malloc(sizeof(double) * width);
    }
  }
  species->delay_val = delay_val;
}

void mySpecies_setLocatingCompartment(mySpecies *species, myCompartment *compartment) {
  species->locating_compartment = compartment;
}

void mySpecies_setDependingRule(mySpecies *species, myRule *rule) {
  species->depending_rule = rule;
}

