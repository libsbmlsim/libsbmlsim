/*
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

#if defined(_MSC_VER) || defined(__STRICT_ANSI__)
#include "libsbmlsim/my_getopt.h"
#else
#include <unistd.h>
#endif

void usage(char *str) {
  printf("Usage : %s [option] filename(SBML file only)\n", str);
  printf(" -t # : specify simulation time (ex. -t 100 )\n");
  printf(" -s # : specify simulation step (ex. -s 100 )\n");
  printf(" -d # : specify simulation delta (ex. -d 0.01 [default:1/4096])\n");
  printf("        dt is calculated in (delta)*(time)/(step)\n");
  printf(" -l   : use lazy method for integration\n");
  printf(" -n   : do not use lazy method\n");
  printf(" -a   : print Species Value in Amount\n");
  printf(" -A # : specify absolute tolerance for variable stepsize (ex. -A 1e-03 [default:1e-09])\n");
  printf(" -R # : specify relative tolerance for variable stepsize (ex. -R 0.1   [default:1e-06])\n");
  printf(" -M # : specify the max change rate of stepsize (ex. -M 1.5 [default:2.0])\n");
  printf(" -B   : use bifurcation analysis \n");
  printf(" -m # : specify numerical integration algorithm (ex. -m 3 )\n");
  printf("        1: Runge-Kutta\n");
  printf("        2: AM1 & BD1 (implicit Euler)\n");
  printf("        3: AM2 (Crank Nicolson)\n");
  printf("        4: AM3\n");
  printf("        5: AM4\n");
  printf("        6: BD2\n");
  printf("        7: BD3\n");
  printf("        8: BD4\n");
  printf("        9: AB1 (explicit Euler)\n");
  printf("       10: AB2\n");
  printf("       11: AB3\n");
  printf("       12: AB4\n");
  printf("       13: Runge-Kutta-Fehlberg\n");
  printf("       14: Cash-Karp\n");
  exit(1);
}

int main(int argc, char *argv[]){
  SBMLDocument_t *d;
  Model_t *m;
  /*  Variables for getopt() */
  int ch;
  extern char *optarg;
  extern int optind;
  unsigned int error_num;

  char *myname;
  boolean use_lazy_method = 2;
  int print_amount = 0;

  double sim_time = 0;
  int step = 0;
  double delta = 1.0/4096;
  double dt = 0;
  int print_interval = 0;

  double atol = ABSOLUTE_ERROR_TOLERANCE;
  double rtol = RELATIVE_ERROR_TOLERANCE;

  char buf1[256], buf2[256], buf3[256];
  char *tmp;
  int method;
  boolean is_explicit;
  char *method_name;
  int method_key = -1;

  /*for bifurcation analysis*/
  boolean use_bifurcation_analysis = false;

  /*for variable step-size numerical integration*/
  boolean use_variable_stepsize = false;
  double facmax = DEFAULT_FACMAX;

  myResult *rtn;

  myname = argv[0];
  while ((ch = getopt(argc, argv, "t:s:d:m:A:R:M:lnaB")) != -1){
    switch (ch) {
      case 't':
        sim_time = atof(optarg);
        break;
      case 's':
        step = atoi(optarg);
        break;
      case 'd':
        delta = atof(optarg);
        break;
      case 'm':
        method_key = atoi(optarg);
        break;
      case 'l':
        use_lazy_method = true;
        break;
      case 'n':
        use_lazy_method = false;
        break;
      case 'a':
        print_amount = 1;
        break;
      case 'A':
		  atol = atof(optarg);
		  break;
      case 'R':
		  rtol = atof(optarg);
		  break;
      case 'B':
		  use_bifurcation_analysis = 1;
		  break;
      case 'M':
		  facmax = atof(optarg);
		  break;
      default:
        usage(myname);
    }
  }
  argc -= optind;
  argv += optind;

  /* get SBML Model */
  if(argc < 1){
    usage(myname);
  }
  d = readSBML(argv[0]);

  error_num = SBMLDocument_getNumErrors(d);
  if (error_num > 0) {
    unsigned int i;
    for (i = 0; i < error_num; i++) {
      const SBMLError_t *err = SBMLDocument_getError(d, i);
      if (XMLError_isError((XMLError_t *)err) || XMLError_isFatal((XMLError_t *)err)) {
        printf("Input file [%s] is not an appropriate SBML file\n", argv[0]);
        exit(1);
      }
    }
  }
  m = SBMLDocument_getModel(d);

  /* determine sim_time */
  if(sim_time == 0){
    while(1){
      printf("Simulation time : ");
      tmp = fgets(buf1, 256, stdin);
      chomp(buf1);
      if(str_is_number(buf1)){
        break;
      }
      printf("not a number!\n");
    }
    sscanf(buf1, "%lf", &sim_time);
  }
  /* determine step */
  if(step == 0){
    while(1){
      printf("simulation step : ");
      tmp = fgets(buf1, 256, stdin);
      chomp(buf1);
      if(str_is_number(buf1)){
        break;
      }
      printf("not a number!\n");
    }
    sscanf(buf1, "%d", &step);
  }

  /* CUI */
  if (method_key == -1) {
    while(1){
      printf("select neumerical integration method\n");
      printf("simulate with\n");
      printf("Runge-Kutta : press \"1\"\n");
      printf("AM1 & BD1 (implicit Euler) : press \"2\"\n");
      printf("AM2 (Crank Nicolson) : press \"3\"\n");
      printf("AM3 : press \"4\"\n");
      printf("AM4 : press \"5\"\n");
      printf("BD2 : press \"6\"\n");
      printf("BD3 : press \"7\"\n");
      printf("BD4 : press \"8\"\n");
      printf("AB1 (explicit Euler) : press \"9\"\n");
      printf("AB2 : press \"10\"\n");
      printf("AB3 : press \"11\"\n");
      printf("AB4 : press \"12\"\n");
	  printf("RKF : press \"13\"\n");
	  printf("CK  : press \"14\"\n");
      tmp = fgets(buf2, 256, stdin);
      method_key = atoi(buf2);
      if (method_key < 1 || method_key > 14) {
        printf("Invalid Input!\nSelect and input the number \"1~12\"");
      } else {
        break;
      }
    }
  }
  switch(method_key) {
    case 1: /*  Runge-Kutta */
      method = MTHD_RUNGE_KUTTA;
      method_name = MTHD_NAME_RUNGE_KUTTA;
      break;
    case 2: /*  Backward-Euler */
      method = MTHD_BACKWARD_EULER;
      method_name = MTHD_NAME_BACKWARD_EULER;
      break;
    case 3: /*  Crank-Nicolson */
      method = MTHD_CRANK_NICOLSON;
      method_name = MTHD_NAME_CRANK_NICOLSON;
      break;
    case 4: /*  Adams-Moulton 3 */
      method = MTHD_ADAMS_MOULTON_3;
      method_name = MTHD_NAME_ADAMS_MOULTON_3;
      break;
    case 5: /*  Adams-Moultion 4 */
      method = MTHD_ADAMS_MOULTON_4;
      method_name = MTHD_NAME_ADAMS_MOULTON_4;
      break;
    case 6: /*  Backward-Difference 2 */
      method = MTHD_BACKWARD_DIFFERENCE_2;
      method_name = MTHD_NAME_BACKWARD_DIFFERENCE_2;
      break;
    case 7: /*  Backward-Difference 3 */
      method = MTHD_BACKWARD_DIFFERENCE_3;
      method_name = MTHD_NAME_BACKWARD_DIFFERENCE_3;
      break;
    case 8: /*  Backward-Difference 4 */
      method = MTHD_BACKWARD_DIFFERENCE_4;
      method_name = MTHD_NAME_BACKWARD_DIFFERENCE_4;
      break;
    case 9: /*  Euler (Adams-Bashforth) */
      method = MTHD_EULER;
      method_name = MTHD_NAME_EULER;
      break;
    case 10: /*  Adams-Bashforth 2 */
      method = MTHD_ADAMS_BASHFORTH_2;
      method_name = MTHD_NAME_ADAMS_BASHFORTH_2;
      break;
    case 11: /*  Adams-Bashforth 3 */
      method = MTHD_ADAMS_BASHFORTH_3;
      method_name = MTHD_NAME_ADAMS_BASHFORTH_3;
      break;
    case 12: /*  Adams-Bashforth 4 */
      method = MTHD_ADAMS_BASHFORTH_4;
      method_name = MTHD_NAME_ADAMS_BASHFORTH_4;
      break;
    case 13: /*  Runge-Kutta-Fehlberg 5 */
      method = MTHD_RUNGE_KUTTA_FEHLBERG_5;
      method_name = MTHD_NAME_RUNGE_KUTTA_FEHLBERG_5;
      use_variable_stepsize = 1;
      break;
    case 14: /*  Cash-Karp */
      method = MTHD_CASH_KARP;
      method_name = MTHD_NAME_CASH_KARP;
      use_variable_stepsize = 1;
      break;
    default:
      method = MTHD_RUNGE_KUTTA;
      method_name = MTHD_NAME_RUNGE_KUTTA;
      break;
  }
  is_explicit = method % 10;

  /* simulation */
  if (!use_variable_stepsize) {
  if (!is_explicit) {
    if (use_lazy_method == 2) {
      while(1){
        printf("use lazy mode?\nlazy mode:using jacobian continuously while the solutions of equations are converging in newton method.\nyes(y) or no(n)\n");
        tmp = fgets(buf3, 256, stdin);
        if(strcmp(buf3, "y\n") == 0 || strcmp(buf3, "yes\n") == 0) {
          use_lazy_method = true;
          break;
        } else if(strcmp(buf3, "n\n") == 0 || strcmp(buf3, "no\n") == 0) {
          use_lazy_method = false;
          break;
        }
        printf("Invalid Input!!!\nSelect yes(y) or no(n)\n");
      }
    }
    if(use_lazy_method == true) {
      printf("  simulate with lazy mode\n");
    }
  }
  }
  /* Run simulation */
  if (!use_variable_stepsize) {
	  /* calculate simulation condition */
	  dt = delta*(sim_time/step);
  }else{
	  /* calculate simulation condition */
	  dt = sim_time/step;
  }
	  print_interval = (int)(1/delta);
	  printf("  time:%g step:%d dt:%f\n", sim_time, step, dt);
	  rtn = simulateSBMLModel(m, sim_time, dt, print_interval, print_amount, method, use_lazy_method, atol, rtol, facmax);

  /* write CSV */
  if (rtn == NULL) {
    printf("Returned result is NULL\n");
  } else {
    write_csv(rtn, "out.csv"); /*  for SBML test suite */
    /* for more generic simulator
       write_separate_result(rtn,
       "./simulation_results/species_result.dat",
       "./simulation_results/parameter_result.dat",
       "./simulation_results/compartment_result.dat");
       */
  }

  /* free */
  SBMLDocument_free(d);
  free_myResult(rtn);
  return 0;
}
