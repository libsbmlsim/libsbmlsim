#include "libsbmlsim/libsbmlsim.h"

int main(int argc, char *argv[]){
  SBMLDocument_t *d;
  Model_t *m;
  myResult result;
  myResult *rtn;

  double sim_time = 4000;
  //double dt = 0.0025 * sim_time / 50;
  double dt = 0.1;
  int print_interval = 100;
  int print_amount = 1;
  int method = 1;
  //int use_lazy_method = 0;
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
  // RK
  rtn = simulateSBMLModel(m, &result, sim_time, dt, print_interval, print_amount, method, use_lazy_method);
  // AM2 (w/o lazy)
  //rtn = simulateSBMLModel(m, &result, sim_time, dt, print_interval, print_amount, 3, 0);
  // AM2 (w/ lazy)
  //rtn = simulateSBMLModel(m, &result, sim_time, dt, print_interval, print_amount, 3, 1);
  if (rtn == NULL) {
    printf("Returned result is NULL\n");
  } else {
    //    print_result(rtn);
    write_result(rtn, "test.dat");
    //    write_csv(rtn, "test.csv");
    write_separate_result(rtn, "./simulation_results/species_result_exp10.dat",
        "./simulation_results/parameter_result_exp10.dat",
        "./simulation_results/compartment_result_exp10.dat");
  }

  SBMLDocument_free(d);
  return 0;
}
