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
#include "libsbmlsim/allocated_memory.h"
#include <stdlib.h>

allocated_memory *allocated_memory_create() {
  allocated_memory *mem = (allocated_memory *)malloc(sizeof(allocated_memory));
  mem->num_of_allocated_memory = 0;
  return mem;
}

void allocated_memory_free(allocated_memory *mem) {
  unsigned int i;

  if (mem == NULL) {
    return;
  }

  for (i = 0; i < mem->num_of_allocated_memory; i++) {
    free(mem->memory[i]);
  }
  free(mem);
}

