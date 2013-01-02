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
#include "libsbmlsim/libsbmlsim.h"

int64_t factorial(int n){
  int i;
  int64_t ans;
  ans = 1;
  for(i=1; i<=n; i++){
    ans *= i;
  }
  return ans;
}

double my_fmax(double a, double b) {
  return a < b ? b : a;
}

double my_fmin(double a, double b) {
  return b < a ? b : a;
}

double my_asinh(double x) {
#if defined(_MSC_VER) || defined(__STRICT_ANSI__)
  if(x == 0.0) {
    return 0.0;
  }
  if (x > 0.0) {
    return log(x + sqrt(x * x + 1));
  } else {
    return -log(-x + sqrt(x * x + 1));
  }
#else
  return asinh(x);
#endif
}

int my_isnan(double x) {
#if defined(_MSC_VER)
  return _isnan(x);
#else
  return isnan(x);
#endif
}
