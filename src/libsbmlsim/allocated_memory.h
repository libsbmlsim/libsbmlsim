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
#ifndef LibSBMLSim_AllocatedMemory_h
#define LibSBMLSim_AllocatedMemory_h

#include "typedefs.h"
#include "common.h"

struct _allocated_memory {
  double *memory[MAX_ALLOCATED_MEMORY];
  unsigned int num_of_allocated_memory;
};

allocated_memory *allocated_memory_create();
void allocated_memory_free(allocated_memory *mem);

#endif /* LibSBMLSim_AllocatedMemory_h */
