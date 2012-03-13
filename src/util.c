#include "libsbmlsim/libsbmlsim.h"

int get_end_cycle(double sim_time, double dt) {
  int r = sim_time / dt;
  return r;
}
