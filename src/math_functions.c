#include <stdio.h> 
#include <math.h>

long long int factorial(int n){
  int i;
  long long int ans;
  ans = 1;
  for(i=1; i<=n; i++){
    ans *= i;
  }
  return ans;
}
