#include "sbml/common/common.h"
#include "sbml/SBMLReader.h"
#include "sbml/SBMLTypes.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include "header.h"

/*AはN*Nの正方行列*/
/*pはピボット交換した行を格納している*/
/*p[0] = 4 ならば、0行目は4行目と交換している*/
/*同時にp[4] = 0 となっている*/
/*b は定数列ベクトル Ax=b*/
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
  /*## Ly=b を解いてyを求める LはAの下三角 ##*/
  /*## 前進代入 ##*/
  for(j=0; j<N; j++){
    sum = b[j];
    for(i=0;i<j;i++){
      sum -= A[j][i]*b[i];
    }
    b[j] = sum;
  }
  /*## Ux=y を解いてx(b)を求める UはAの上三角 ##*/
  /*## 後退代入 ##*/
  for(j=N-1; j>=0; j--){
    sum = b[j];
    for(i=N-1; i>j; i--){
      sum -= A[j][i]*b[i];
    }
    b[j] = sum / A[j][j];
  }
  
  return 0;
}
