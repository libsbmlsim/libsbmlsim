/**
 * <!--------------------------------------------------------------------------
 * This file is part of libSBMLSim.  Please visit
 * http://fun.bio.keio.ac.jp/software/libsbmlsim/ for more
 * information about libSBMLSim and its latest version.
 *
 * Copyright (C) 2011-2016 by the Keio University, Yokohama, Japan
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution.
 * ---------------------------------------------------------------------- -->*/
#ifndef LibSBMLSim_MyEventAssignment_h
#define LibSBMLSim_MyEventAssignment_h

#include "typedefs.h"
#include "common.h"
#include "equation.h"
#include "mySpecies.h"
#include "myParameter.h"
#include "myCompartment.h"
#include "mySpeciesReference.h"
#include <sbml/SBMLTypes.h>

struct _myEventAssignment {
  EventAssignment_t *origin;
  equation *eq; /* assignment equation */
  mySpecies *target_species;
  myParameter *target_parameter;
  myCompartment *target_compartment;
  mySpeciesReference *target_species_reference;
};

myEventAssignment *myEventAssignment_create();
void myEventAssignment_initWithEvent(myEventAssignment *assign, Event_t *event, int index);
void myEventAssignment_free(myEventAssignment *assign);

void myEventAssignment_initTarget(
    myEventAssignment *assign,
    mySpecies **species, unsigned int num_of_species,
    myParameter **parameters, unsigned int num_of_parameters,
    myCompartment **compartments, unsigned int num_of_compartments,
    myReaction **reactions, unsigned int num_of_reactions);

EventAssignment_t *myEventAssignment_getOrigin(myEventAssignment *assign);
void myEventAssignment_setTargetSpecies(myEventAssignment *assign, mySpecies *species);
void myEventAssignment_setTargetParameter(myEventAssignment *assign, myParameter *parameter);
void myEventAssignment_setTargetCompartment(myEventAssignment *assign, myCompartment *compartment);
void myEventAssignment_setTargetSpeciesReference(myEventAssignment *assign, mySpeciesReference *ref);

#endif /* LibSBMLSim_MyEventAssignment_h */
