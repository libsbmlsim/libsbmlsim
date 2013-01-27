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
#include "libsbmlsim/myDelay.h"
#include <stdlib.h>
#include <sbml/SBMLTypes.h>

myDelay *myDelay_create() {
  myDelay *delay = (myDelay *)malloc(sizeof(myDelay));
  delay->origin = NULL;
  delay->eq = NULL;
  return delay;
}

void myDelay_initWithOrigin(myDelay *delay, Delay_t *origin) {
  delay->origin = origin;
  delay->eq = equation_create();
}

void myDelay_free(myDelay *delay) {
  if (delay == NULL) {
    return;
  }

  if (delay->eq != NULL) {
    equation_free(delay->eq);
  }
  free(delay);
}

Delay_t *myDelay_getOrigin(myDelay *delay) {
  return delay->origin;
}

equation *myDelay_getEquation(myDelay *delay) {
  return delay->eq;
}

