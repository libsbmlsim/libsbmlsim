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
#ifndef LibSBMLSim_Variable_h
#define LibSBMLSim_Variable_h

#include "typedefs.h"

struct _variable {
  char *id;
  double value;
};

variable *variable_create(const char *id, double value);
void variable_free(variable *var);
const char *variable_getId(variable *var);
void variable_setValue(variable *var, double value);
double variable_getValue(variable *var);

#endif /* LibSBMLSim_Variable_h */
