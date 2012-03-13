#include "libsbmlsim/libsbmlsim.h"

/* LU decomposition with pivoting
 * A = N*N (matrix)
 * p contains an index of rows where pivot is exchanged.
 * (ex.) p[0] = 4 ... row 0 and 4 are exchanged.
 *                    at the same time, p[4] = 0.
 */
int lu_decomposition(double **A, int *p, int N){
  float Eps;
  int i, j, t;
  int mj;
  float tmp;
  int   tmp2;
  Eps = 1e-10;

  for(t=0; t<N; t++){
    /* PIVOT */
    mj = t;
    for(j=t+1;j<N;j++){
      if(fabs(A[j][t]) > fabs(A[mj][t])){
        mj = j;
      }
    }
    if(mj != t){
      for(i=0; i<N; i++){
        tmp = A[mj][i];
        A[mj][i] = A[t][i];
        A[t][i]  = tmp;
      }
      tmp2 = p[mj];
      p[mj] = p[t];
      p[t] = tmp2;
    }
    /* END PIVOT */
    if(fabs(A[t][t]) < Eps){
      dbg_printf("A is singular matrix \n");
      return 0;
    }
    for(j=t+1; j<N; j++){
      A[j][t] = A[j][t]/A[t][t]; 
      for(i=t+1; i<N; i++){
        A[j][i] = A[j][i] - A[j][t]*A[t][i];
      }
    }
  }
  return 1;
}
