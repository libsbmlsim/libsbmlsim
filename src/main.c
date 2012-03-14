#include <unistd.h>
#include "libsbmlsim/libsbmlsim.h"

void usage(char *str) {
  printf("Usage : %s [option] filename(SBML file only)\n", str);
  printf(" -t # : specify simulation time (ex. -t 100 )\n");
  printf(" -s # : specify simulation step (ex. -s 100 )\n");
  printf(" -d # : specify simulation delta (ex. -d 0.01 [default:1/4092])");
  printf("        dt is calculated in (delta)*(time)/(step)\n");
  printf(" -l   : use lazy method for integration\n");
  printf(" -n   : do not use lazy method\n");
  printf(" -a   : print Species Value in Amount\n");
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
  exit(1);
}

/* Extended SBML_simulator_1028  */
/* improved in main.c optget */
/* processing for compartment is added */
int main(int argc, char *argv[]){
  SBMLDocument_t *d;
  Model_t *m;

  /*  Valuables for getopt() */
  int ch;
  extern char *optarg;
  extern int optind, opterr;

  char *myname;
  boolean use_lazy_method = -1;
  int print_amount = 0;

  double sim_time = 0;
  int step = 0;
  double delta = 1.0/4092;
  double dt = 0;
  int print_interval = 0;

  char buf1[256], buf2[256], buf3[256];
  int method;
  int order = 0;
  int is_explicit;
  char *method_name;
  int method_key = -1;

  allocated_memory *mem;
  mem = (allocated_memory*)malloc(sizeof(allocated_memory));
  mem->num_of_allocated_memory = 0;

  copied_AST *cp_AST;
  cp_AST = (copied_AST*)malloc(sizeof(copied_AST));
  cp_AST->num_of_copied_AST = 0;

  myname = argv[0];
  while ((ch = getopt(argc, argv, "t:s:d:m:lna")) != -1){
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
  /*   if(SBMLDocument_getNumErrors(d) > 0){ */
  /*     printf("Input file [%s] is not an appropriate SBML file\n", argv[0]); */
  /*     exit(1); */
  /*   } */
  m = SBMLDocument_getModel(d);

  /* determine sim_time */
  if(sim_time == 0){
    while(1){
      printf("Simulation time : ");
      fgets(buf1, 256, stdin);
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
      fgets(buf1, 256, stdin);
      chomp(buf1);
      if(str_is_number(buf1)){
        break;
      }
      printf("not a number!\n");
    }
    sscanf(buf1, "%d", &step);
  }

  /* calculate simulation condition */
  dt = delta*(sim_time/step);
  print_interval = (int)(1/delta);
  printf("  time:%g step:%d dt:%lf\n", sim_time, step, dt);

  /* time in simulation */
  double time = 0;
  /* prepare mySpecies */
  int num_of_species = Model_getNumSpecies(m);
  mySpecies *mySp[num_of_species];
  /* prepare myParameters */
  int num_of_parameters = Model_getNumParameters(m);
  myParameter *myParam[num_of_parameters];
  /* prepare myCompartments */
  int num_of_compartments = Model_getNumCompartments(m);
  myCompartment *myComp[num_of_compartments];
  /* prepare myReactions */
  int num_of_reactions = Model_getNumReactions(m);
  myReaction *myRe[num_of_reactions];
  /* prepare myRules */
  int num_of_rules = Model_getNumRules(m);
  myRule *myRu[num_of_rules];
  /* prepare myEvents */
  int num_of_events = Model_getNumEvents(m);  
  myEvent *myEv[num_of_events];
  /* prepare myInitial Assignments */
  int num_of_initialAssignments = Model_getNumInitialAssignments(m);
  myInitialAssignment *myInitAssign[num_of_initialAssignments];
  /* prepare myAlgebraicEquations */
  myAlgebraicEquations *myAlgEq = NULL;
  /* prepare timeVariantAssignments */
  timeVariantAssignments *timeVarAssign = NULL;
  /* create myObjects */
  create_mySBML_objects(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, &myAlgEq, &timeVarAssign, sim_time, dt, &time, mem, cp_AST);
  /* prepare myResult */
  myResult result;
  /* create myResult */
  create_myResult_content(m, &result, mySp, myParam, myComp, sim_time, dt, print_interval);
  /* prepare return value */
  myResult* rtn;
  if(myAlgEq == NULL){
    dbg_printf("myAlgEq is NULL\n");
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

      fgets(buf2, 256, stdin);
      method_key = atoi(buf2);
      if (method_key < 1 || method_key > 12) {
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
    case 6: /*  Backward-Differentiation 2 */
      method = MTHD_BACKWARD_DIFFERENTIATION_2;
      method_name = MTHD_NAME_BACKWARD_DIFFERENTIATION_2;
      break;
    case 7: /*  Backward-Differentiation 3 */
      method = MTHD_BACKWARD_DIFFERENTIATION_3;
      method_name = MTHD_NAME_BACKWARD_DIFFERENTIATION_3;
      break;
    case 8: /*  Backward-Differentiation 4 */
      method = MTHD_BACKWARD_DIFFERENTIATION_4;
      method_name = MTHD_NAME_BACKWARD_DIFFERENTIATION_4;
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
    default:
      method = MTHD_RUNGE_KUTTA;
      method_name = MTHD_NAME_RUNGE_KUTTA;
      break;
  }
  order = method / 10;
  is_explicit = method % 10;
  dbg_printf("simulate with %s\n", method_name);

  /* simulation */
  if (is_explicit == 1) {
    rtn = simulate_explicit(m, &result, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval, &time, order, print_amount, mem);
  }else{
    if (use_lazy_method == -1) {
      while(1){
        printf("use lazy mode?\nlazy mode:using jacobian continuously while the solutions of equations are converging in newton method.\nyes(y) or no(n)\n");
        fgets(buf3, 256, stdin);
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
    rtn = simulate_implicit(m, &result, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval, &time, order, use_lazy_method, print_amount, mem);
  }

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

  /* print result list */
  /* print_result_list(m, mySp, myParam, myComp); */

  /* free */
  printf("  free all allocated memory\n");
  free_mySBML_objects(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, mem, cp_AST);
  SBMLDocument_free(d);
  return 0;
}
