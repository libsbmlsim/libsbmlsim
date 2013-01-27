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
#ifndef LibSBMLSim_MyEventAssignment_h
#define LibSBMLSim_MyEventAssignment_h

#include "typedefs.h"
#include <sbml/SBMLTypes.h>

struct _myEventAssignment {
  EventAssignment_t *origin;
  equation *eq; /* assignment equation */
  mySpecies *target_species;
  myParameter *target_parameter;
  myCompartment *target_compartment;
  mySpeciesReference *target_species_reference;
};

#endif /* LibSBMLSim_MyEventAssignment_h */
