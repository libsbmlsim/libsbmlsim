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
#ifndef LibSBMLSim_MyDelay_h
#define LibSBMLSim_MyDelay_h

#include "typedefs.h"
#include "boolean.h"
#include "equation.h"
#include <sbml/SBMLTypes.h>

struct _myDelay {
  Delay_t *origin;
  equation *eq;
};

myDelay *myDelay_create();
void myDelay_initWithOrigin(myDelay *delay, Delay_t *origin);
void myDelay_free(myDelay *delay);

Delay_t *myDelay_getOrigin(myDelay *delay);
equation *myDelay_getEquation(myDelay *delay);

#endif /* LibSBMLSim_MyDelay_h */
