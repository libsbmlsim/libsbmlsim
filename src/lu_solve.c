#include "libsbmlsim/libsbmlsim.h"

/* 
 * A = N*N (matrix).
 * p contains an index of rows where pivot is exchanged.
 * b is vector of contant column: Ax = b.
 * (ex.) p[0] = 4 ... row 0 and 4 are exchanged.
 *                    at the same time, p[4] = 0.
 */
int lu_solve(double **A, int *p, int N, double *b){
  float sum;
  int i, j;
  float tmp;

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
