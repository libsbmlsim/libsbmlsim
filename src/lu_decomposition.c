#include "sbml/common/common.h"
#include "sbml/SBMLReader.h"
#include "sbml/SBMLTypes.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include "header.h"

/*ピボット交換ありのLU分解*/
/*AはN*Nの正方行列*/
/*pはピボット交換した行を格納している*/
/*p[0] = 4 ならば、0行目は4行目と交換している*/
/*同時にp[4] = 0 となっている*/
int lu_decomposition(double **A, int *p, int N){
  float Eps;
  int i, j, t;
  int mj;
  float tmp;
  int   tmp2;
  Eps = 1e-10;

  for(t=0; t<N; t++){
    /*## PIVOT ##*/
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
    /*## END PIVOT ##*/
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
