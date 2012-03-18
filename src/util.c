#include <time.h>
#include "libsbmlsim/libsbmlsim.h"

int get_end_cycle(double sim_time, double dt) {
  int r = (int)(sim_time / dt);
  return r;
}

void set_seed(void){
  srand((unsigned)time(NULL));
}

char* dupstr(const char *str)
{
  char *copy = NULL;
  if (str) {
    copy = malloc(strlen(str)+1);
    if (copy) {
#ifdef _MSC_VER
      strcpy_s(copy, strlen(str)+1, str);
#else
      strcpy(copy, str);
#endif
    }
  }
  return copy;
}
