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
#include<stdio.h>
#include<stdlib.h>

/* Total bytes allocated */
extern int total_allocated;

/* Memory alignment is important */
typedef union { double d; struct {size_t n; char *file; int line;} s; } Site;

void* debug_malloc(size_t n, char *file, int line) {
  char *rp;
  rp = (char*)malloc(sizeof(Site)+n);
  total_allocated += n;
  ((Site*)rp)->s.n = n;
  ((Site*)rp)->s.file = file;
  ((Site*)rp)->s.line = line;
  return (void*)(rp + sizeof(Site));
}

void* debug_calloc(size_t c, size_t n, char *file, int line) {
  char *rp;
  rp = (char*)calloc(sizeof(Site)+n*c, sizeof(char));
  total_allocated += n*c;
  ((Site*)rp)->s.n = n*c;
  ((Site*)rp)->s.file = file;
  ((Site*)rp)->s.line = line;
  return (void*)(rp + sizeof(Site));
}

void* debug_realloc(void *ptr, size_t n, char *file, int line) {
  char *rp;
  char *tmp;
  tmp = ((char*)ptr) - sizeof(Site);
  total_allocated -= ((Site*)tmp)->s.n;
  rp = (char*)realloc(tmp, sizeof(Site)+n);
  total_allocated += n;
  ((Site*)rp)->s.n = n;
  ((Site*)rp)->s.file = file;
  ((Site*)rp)->s.line = line;
  return (void*)(rp + sizeof(Site));
}

void debug_free(void *p, char *file, int line) {
  char *rp;
  rp = ((char*)p) - sizeof(Site);
  total_allocated -= ((Site*)rp)->s.n;
  ((Site*)rp)->s.file = file;
  ((Site*)rp)->s.line = line;
  free(rp);
}
