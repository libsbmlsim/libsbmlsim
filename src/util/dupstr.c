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
#include <stdlib.h>
#include <string.h>

char* dupstr(const char *str)
{
  char *copy = NULL;
  if (str) {
    copy = malloc(strlen(str)+1);
    if (copy) {
#ifdef _MSC_VER
      strcpy_s(copy, strlen(str)+1, str);
#else
      strcpy(copy, str);
#endif
    }
  }
  return copy;
}
