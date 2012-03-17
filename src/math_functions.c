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
