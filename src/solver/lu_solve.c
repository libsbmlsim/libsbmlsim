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
#include "../libsbmlsim/libsbmlsim.h"

/* 
 * A = N*N (matrix).
 * p contains an index of rows where pivot is exchanged.
 * b is vector of contant column: Ax = b.
 * (ex.) p[0] = 4 ... row 0 and 4 are exchanged.
 *                    at the same time, p[4] = 0.
 */
int lu_solve(double **A, int *p, int N, double *b){
  double sum;
  int i, j;
  double tmp;

  for(j=0; j<N; j++){
    if(j != p[j] && j < p[j]){
      tmp = b[j];
      b[j] = b[p[j]];
      b[p[j]]= tmp;
    }
  }
  /* Solve Ly = b, and obtain y. L: lower triangular matrix of A */
  /* Forward substitution */
  for(j=0; j<N; j++){
    sum = b[j];
    for(i=0;i<j;i++){
      sum -= A[j][i]*b[i];
    }
    b[j] = sum;
  }
  /* Solve Ux = y, and obtain x(b). U: upper triangular matrix of A */
  /* Backward substitution */
  for(j=N-1; j>=0; j--){
    sum = b[j];
    for(i=N-1; i>j; i--){
      sum -= A[j][i]*b[i];
    }
    b[j] = sum / A[j][j];
  }

  return 0;
}
