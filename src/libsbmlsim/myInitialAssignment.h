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
#ifndef LibSBMLSim_MyInitialAssignment_h
#define LibSBMLSim_MyInitialAssignment_h

#include "typedefs.h"
#include "boolean.h"
#include "equation.h"
#include "mySpecies.h"
#include "myParameter.h"
#include "myCompartment.h"
#include "mySpeciesReference.h"
#include <sbml/SBMLTypes.h>

struct _myInitialAssignment {
  InitialAssignment_t *origin;
  equation *eq;
  mySpecies *target_species;
  myParameter *target_parameter;
  myCompartment *target_compartment;
  mySpeciesReference *target_species_reference;
};

myInitialAssignment *myInitialAssignment_create();
void myInitialAssignment_initWithModel(myInitialAssignment *assign, Model_t *model, int index);
void myInitialAssignment_free(myInitialAssignment *assign);
void myInitialAssignment_initTarget(
    myInitialAssignment *assign,
    mySpecies **species, unsigned int num_of_species,
    myParameter **parameters, unsigned int num_of_parameters,
    myCompartment **compartments, unsigned int num_of_compartments);
void myInitialAssignment_initTargetSpeciesReference(
    myInitialAssignment *assign,
    myReaction **reactions, unsigned int num_of_reactions);

InitialAssignment_t *myInitialAssignment_getOrigin(myInitialAssignment *assign);
void myInitialAssignment_setTargetSpecies(myInitialAssignment *assign, mySpecies *species);
void myInitialAssignment_setTargetParameter(myInitialAssignment *assign, myParameter *parameter);
void myInitialAssignment_setTargetCompartment(myInitialAssignment *assign, myCompartment *compartment);
void myInitialAssignment_setTargetSpeciesReference(myInitialAssignment *assign, mySpeciesReference *ref);

#endif /* LibSBMLSim_MyInitialAssignment_h */
