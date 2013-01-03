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
#include <ctype.h>
#include "../libsbmlsim/libsbmlsim.h"

boolean str_is_number(const char *str)
{
  boolean isfloat = false;
  if (*str == '\0')
    return false;
  while (*str != '\0') {
    if (isdigit(*str)) {
      ;
    } else if (*str == '.') {
      if (isfloat)
        return false;
      isfloat = true;
    } else {
      return false;
    }
    str++;
  }
  return true;
}

