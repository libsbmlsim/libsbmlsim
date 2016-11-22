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
#include "typedefs.h"
#include "libsbmlsim/libsbmlsim.h"
#include <stdlib.h>
#include <string.h>

variable *variable_create(const char *id, double value) {
  variable *var = (variable *)malloc(sizeof(variable));
  int idlen = strlen(id);
  var->id = (char *)malloc(sizeof(char) * (idlen + 1));
  strncpy(var->id, id, idlen);
  var->value = value;
  return var;
}

void variable_free(variable *var) {
  free(var->id);
  free(var);
}

const char *variable_getId(variable *var) {
  return var->id;
}

void variable_setValue(variable *var, double value) {
  var->value = value;
}

double variable_getValue(variable *var) {
  return var->value;
}

