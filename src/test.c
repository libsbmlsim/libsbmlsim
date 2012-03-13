#include "libsbmlsim/libsbmlsim.h"

int main(int argc, char *argv[]) {
  SBMLDocument_t *d;
  Model_t *m;
  myResult result;
  myResult *rtn;

  double sim_time = 4000;
  double dt = 0.1;
  int print_interval = 100;
  int print_amount = 1;
  int method = MTHD_RUNGE_KUTTA;
  boolean use_lazy_method = false;

  if (argc < 2) {
    printf("Input SBML file is not specified.\n  Usage: %s sbml.xml\n", argv[0]);
    exit(1);
  }
  d = readSBML(argv[1]);
  if (SBMLDocument_getNumErrors(d) > 0) {
    printf("Input file [%s] is not an appropriate SBML file\n", argv[1]);
    exit(1);
  }
  m = SBMLDocument_getModel(d);
  rtn = simulateSBMLModel(m, &result, sim_time, dt, print_interval, print_amount, method, use_lazy_method);
  if (rtn == NULL) {
    printf("Returned result is NULL\n");
  } else {
    // print_result(rtn);
    // write_csv(rtn, "test.csv");
    write_result(rtn, "test.dat");
  }

  SBMLDocument_free(d);
  return 0;
}
