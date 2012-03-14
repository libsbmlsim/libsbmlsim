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
