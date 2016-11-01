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
#include "libsbmlsim/myCompartment.h"
#include <stdlib.h>
#include <sbml/SBMLTypes.h>

static unsigned int get_num_of_including_species(Compartment_t *compartment, Model_t *model);

myCompartment *myCompartment_create() {
  myCompartment *compartment = (myCompartment *)malloc(sizeof(myCompartment));
  compartment->origin = NULL;
  compartment->value = 1.0;
  compartment->temp_value = 1.0;
  compartment->k[0] = 0;
  compartment->k[1] = 0;
  compartment->k[2] = 0;
  compartment->k[3] = 0;
  compartment->k[4] = 0;
  compartment->k[5] = 0;
  compartment->delay_val = NULL;
  compartment->delay_val_width = 0;
  compartment->delay_val_length = 0;
  compartment->depending_rule = NULL;
  compartment->prev_val[0] = 1.0;
  compartment->prev_val[1] = 1.0;
  compartment->prev_val[2] = 1.0;
  compartment->prev_k[0] = 0;
  compartment->prev_k[1] = 0;
  compartment->prev_k[2] = 0;
  compartment->including_species = NULL;
  compartment->num_of_including_species = 0;
  return compartment;
}

void myCompartment_initWithModel(myCompartment *compartment, Model_t *model, int index) {
  Compartment_t *origin = (Compartment_t *)ListOf_get(Model_getListOfCompartments(model), index);
  unsigned int num_of_including_species;

  compartment->origin = origin;
  if (Compartment_isSetSize(origin)) {
    compartment->value = Compartment_getSize(origin);
  } else {
    compartment->value = 1.0;
  }
  compartment->temp_value = compartment->value;
  compartment->prev_val[0] = compartment->value;
  compartment->prev_val[1] = compartment->value;
  compartment->prev_val[2] = compartment->value;

  num_of_including_species = get_num_of_including_species(origin, model);
  compartment->including_species = (mySpecies **)malloc(num_of_including_species * sizeof(mySpecies *));
}

void myCompartment_initDelayVal(myCompartment *compartment, unsigned int length, unsigned int width) {
  unsigned int i;
  compartment->delay_val_width = width;
  compartment->delay_val_length = length;
  compartment->delay_val = (double **)malloc(sizeof(double *) * length);
  for (i = 0; i < length; i++) {
    compartment->delay_val[i] = (double *)malloc(sizeof(double) * width);
  }
}

void myCompartment_free(myCompartment *compartment) {
  unsigned int i;
  if (compartment->delay_val != NULL) {
    for (i = 0; i < compartment->delay_val_length; i++) {
      free(compartment->delay_val[i]);
    }
    free(compartment->delay_val);
  }
  if (compartment->including_species != NULL) {
    free(compartment->including_species);
  }
  free(compartment);
}

void myCompartment_reallocDelayVal(myCompartment *compartment, unsigned int length, unsigned int width) {
  unsigned int i;
  unsigned int old_width = compartment->delay_val_width;
  unsigned int old_length = compartment->delay_val_length;
  double **delay_val = (double **)realloc(compartment->delay_val, sizeof(double *) * length);

  compartment->delay_val_width = width;
  compartment->delay_val_length = length;
  for (i = 0; i < length; i++) {
    if (i < old_length && width != old_width) {
      delay_val[i] = (double *)realloc(compartment->delay_val[i], sizeof(double) * width);
    } else if (i >= old_length) {
      delay_val[i] = (double *)malloc(sizeof(double) * width);
    }
  }
  compartment->delay_val = delay_val;
}

Compartment_t *myCompartment_getOrigin(myCompartment *compartment) {
  return compartment->origin;
}

void myCompartment_setDependingRule(myCompartment *compartment, myRule *rule) {
  compartment->depending_rule = rule;
}

void myCompartment_addIncludingSpecies(myCompartment *compartment, mySpecies *species) {
  unsigned int num = compartment->num_of_including_species;
  compartment->including_species[num] = species;
  compartment->num_of_including_species++;
}

static unsigned int get_num_of_including_species(Compartment_t *compartment, Model_t *model) {
  unsigned int ret = 0;
  const char *cid = Compartment_getId(compartment);
  unsigned int num_of_species = Model_getNumSpecies(model);
  unsigned int i;
  for (i = 0; i < num_of_species; i++) {
    Species_t *s = Model_getSpecies(model, i);
    if (Species_isSetCompartment(s) &&
        strcmp(cid, Species_getCompartment(s)) == 0) {
      ret++;
    }
  }
  return ret;
}

