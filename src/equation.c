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
#include "libsbmlsim/equation.h"
#include <stdlib.h>
#include <sbml/SBMLTypes.h>

equation *equation_create() {
  equation *ret = (equation *)malloc(sizeof(equation));
  ret->math_length = 0;
  return ret;
}

void equation_free(equation *eq) {
  if (eq == NULL) {
    return;
  }
  free(eq);
}
