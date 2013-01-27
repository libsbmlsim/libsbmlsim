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
#ifndef LibSBMLSim_MyRule_h
#define LibSBMLSim_MyRule_h

#include "typedefs.h"
#include "boolean.h"
#include "equation.h"
#include "mySpecies.h"
#include "myParameter.h"
#include "myCompartment.h"
#include "mySpeciesReference.h"
#include <sbml/SBMLTypes.h>

struct _myRule {
  Rule_t *origin;
  equation *eq;
  mySpecies *target_species;
  myParameter *target_parameter;
  myCompartment *target_compartment;
  mySpeciesReference *target_species_reference;
  boolean is_rate;
  boolean is_assignment;
  boolean is_algebraic;
};

myRule *myRule_create();
void myRule_initWithModel(myRule *rule, Model_t *model, int index);
void myRule_free(myRule *rule);

void myRule_initTarget(
    myRule *rule,
    mySpecies **species, unsigned int num_of_species,
    myParameter **parameters, unsigned int num_of_parameters,
    myCompartment **compartments, unsigned int num_of_compartments,
    myReaction **reactions, unsigned int num_of_reactions);

Rule_t *myRule_getOrigin(myRule *rule);
void myRule_setTargetSpecies(myRule *rule, mySpecies *species);
void myRule_setTargetParameter(myRule *rule, myParameter *parameter);
void myRule_setTargetCompartment(myRule *rule, myCompartment *compartment);
void myRule_setTargetSpeciesReference(myRule *rule, mySpeciesReference *ref);
equation *myRule_getEquation(myRule *rule);

#endif /* LibSBMLSim_MyRule_h */
