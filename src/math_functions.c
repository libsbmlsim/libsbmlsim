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

double _asinh(double x) {
#if defined(WIN32) || (__STRICT_ANSI__)
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
