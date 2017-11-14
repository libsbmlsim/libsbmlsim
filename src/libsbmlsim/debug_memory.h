/**
 * <!--------------------------------------------------------------------------
 * This file is part of libSBMLSim.  Please visit
 * http://fun.bio.keio.ac.jp/software/libsbmlsim/ for more
 * information about libSBMLSim and its latest version.
 *
 * Copyright (C) 2011-2017 by the Keio University, Yokohama, Japan
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution.
 * ---------------------------------------------------------------------- -->*/
#ifndef LibSBMLSim_Debug_Memory_h
#define LibSBMLSim_Debug_Memory_h

/* Total bytes allocated */
int total_allocated;

/* wrapper for malloc, free */
void* debug_malloc(size_t, char*, int);
void* debug_calloc(size_t, size_t, char*, int);
void* debug_realloc(void*, size_t, char*, int);
void debug_free(void*, char*, int);
char* debug_strdup(const char*, char*, int);
char* debug_strndup(const char*, size_t, char*, int);

#endif /* LibSBMLSim_Debug_Memory_h */
