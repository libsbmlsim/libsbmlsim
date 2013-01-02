/**
 * <!--------------------------------------------------------------------------
 * This file is part of libSBMLSim.  Please visit
 * http://fun.bio.keio.ac.jp/software/libsbmlsim/ for more
 * information about libSBMLSim and its latest version.
 *
 * Copyright (C) 2011-2013 by the Keio University, Yokohama, Japan
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution.
 * ---------------------------------------------------------------------- -->*/
#include "libsbmlsim/libsbmlsim.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {
  myResult *rtn;

  double sim_time = 25;
  double dt = 0.01;
  int print_interval = 10;
  int print_amount = 1;
  int method = MTHD_RUNGE_KUTTA;
  boolean use_lazy_method = false;

  if (argc < 2) {
    cout << "Input SBML file is not specified.\n  Usage: " << argv[0] << "sbml.xml\n" << endl;
    exit(1);
  }
  rtn = simulateSBMLFromFile(argv[1], sim_time, dt, print_interval, print_amount, method, use_lazy_method);
  if (rtn == NULL) {
    cout << "Returned result is NULL" << endl;
  } else if (myResult_isError(rtn)) {
    cout << "ERROR: " <<  myResult_getErrorMessage(rtn) << endl;
  } else {
    /*  print_result(rtn); */
    /* write_csv(rtn, "cresult.csv"); */
    write_result(rtn, "test.dat");
    cout << "Simulation result is written to test.dat." << endl;
  }
  if (rtn != NULL)
    free_myResult(rtn);
  return 0;
}
