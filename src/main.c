#include "sbml/common/common.h"
#include "sbml/SBMLReader.h"
#include "sbml/SBMLTypes.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include <unistd.h>
#include "header.h"

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
  printf("        1: runge kutta\n");
  printf("        2: AM1 & BD1 (implicit eular)\n");
  printf("        3: AM2 (crank nicolson)\n");
  printf("        4: AM3\n");
  printf("        5: AM4\n");
  printf("        6: BD2\n");
  printf("        7: BD3\n");
  printf("        8: BD4\n");
  printf("        9: AB1 (explicit eular)\n");
  printf("       10: AB2\n");
  printf("       11: AB3\n");
  printf("       12: AB4\n");
  exit(1);
}

//Extended SBML_simulator_1028 
//improved in main.c optget
//processing for compartment is added
int main(int argc, char *argv[]){
  SBMLDocument_t *d;
  Model_t *m;

  // Valuables for getopt()
  int ch;
  extern char *optarg;
  extern int optind, opterr;

  char *myname;
  int use_lazy_method = -1;
  int is_explicit = 0;
  int print_amount = 0;

  double sim_time = 0;
  int step = 0;
  double delta = 1.0/4092;
  double dt = 0;
  int print_interval = 0;

  char buf1[256], buf2[256], buf3[256];
  int order = 0;
  int method = -1;

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
      method = atoi(optarg);
      break;
    case 'l':
      use_lazy_method = 1;
      break;
    case 'n':
      use_lazy_method = 0;
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

  //get SBML Model
  if(argc < 1){
    usage(myname);
  }
  d = readSBML(argv[0]);
/*   if(SBMLDocument_getNumErrors(d) > 0){ */
/*     printf("Input file [%s] is not an appropriate SBML file\n", argv[0]); */
/*     exit(1); */
/*   } */
  m = SBMLDocument_getModel(d);
  
  //determine sim_time
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
  //determine step
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

  //calculate simulation condition
  dt = delta*(sim_time/step);
  print_interval = (int)(1/delta);
  printf("  time:%g step:%d dt:%lf\n", sim_time, step, dt);

  //time in simulation
  double time = 0;
  //prepare mySpecies
  int num_of_species = Model_getNumSpecies(m);
  mySpecies *mySp[num_of_species];
  //prepare myParameters
  int num_of_parameters = Model_getNumParameters(m);
  myParameter *myParam[num_of_parameters];
  //prepare myCompartments
  int num_of_compartments = Model_getNumCompartments(m);
  myCompartment *myComp[num_of_compartments];
  //prepare myReactions
  int num_of_reactions = Model_getNumReactions(m);
  myReaction *myRe[num_of_reactions];
  //prepare myRules
  int num_of_rules = Model_getNumRules(m);
  myRule *myRu[num_of_rules];
  //prepare myEvents
  int num_of_events = Model_getNumEvents(m);  
  myEvent *myEv[num_of_events];
  //prepare myInitial Assignments
  int num_of_initialAssignments = Model_getNumInitialAssignments(m);
  myInitialAssignment *myInitAssign[num_of_initialAssignments];
  //prepare myAlgebraicEquations
  myAlgebraicEquations *myAlgEq = NULL;
  //prepare timeVariantAssignments
  timeVariantAssignments *timeVarAssign = NULL;
  //create myObjects
  create_mySBML_objects(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, &myAlgEq, &timeVarAssign, sim_time, dt, &time, mem, cp_AST);
  //prepare myResult
  myResult result;
  //create myResult
  create_myResult_content(m, &result, mySp, myParam, myComp, sim_time, dt, print_interval);
  //prepare return value
  myResult* rtn;
  //
  if(myAlgEq == NULL){
    dbg_printf("myAlgEq is NULL\n");
  }
  //CUI
  if (method == -1) {
    while(1){
      printf("select neumerical integration method\n");
      printf("simulate with\n");
      printf("runge kutta : push \"1\"\n");
      printf("AM1 & BD1 (implicit eular) : push \"2\"\n");
      printf("AM2 (crank nicolson) : push \"3\"\n");
      printf("AM3 : push \"4\"\n");
      printf("AM4 : push \"5\"\n");
      printf("BD2 : push \"6\"\n");
      printf("BD3 : push \"7\"\n");
      printf("BD4 : push \"8\"\n");
      printf("AB1 (explicit eular) : push \"9\"\n");
      printf("AB2 : push \"10\"\n");
      printf("AB3 : push \"11\"\n");
      printf("AB4 : push \"12\"\n");

      fgets(buf2, 256, stdin);
      method = atoi(buf2);
      if (method < 1 || method > 12) {
	printf("Invalid Input!\nSelect and input the number \"1~12\"");
      } else {
	break;
      }
    }
  }
  switch(method) {
  case 1:
    printf("  simulate with runge kutta\n");
    order = 4;
    is_explicit = 1;
    break;
  case 2:
    printf("  simulate with Backward-Eular\n");
    order = 0;
    break;
  case 3:
    printf("  simulate with AM2(Crank-Nicolson)\n");
    order = 1;
    break;
  case 4:
    printf("  simulate with AM3\n");
    order = 2;
    break;
  case 5:
    printf("  simulate with AM4\n");
    order = 3;
    break;
  case 6:
    printf("  simulate with BD2\n");
    order = 4;
    break;
  case 7:
    printf("  simulate with BD3\n");
    order = 5;
    break;
  case 8:
    printf("  simulate with BD4\n");
    order = 6;
    break;
  case 9:
    printf("  simulate with AB1(eular)\n");
    order = 0;
    is_explicit = 1;
    break;
  case 10:
    printf("  simulate with AB2\n");
    order = 1;
    is_explicit = 1;
    break;
  case 11:
    printf("  simulate with AB3\n");
    order = 2;
    is_explicit = 1;
    break;
  case 12:
    printf("  simulate with AB4\n");
    order = 3;
    is_explicit = 1;
    break;
  default:
    printf("  simulate with runge kutta\n");
    order = 4;
    break;
  }

  //simulation
  if (is_explicit == 1) {
    rtn = simulate_explicit(m, &result, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval, &time, order, print_amount, mem);
  }else{
    if (use_lazy_method == -1) {
      while(1){
	printf("use lazy mode?\nlazy mode:using jacobian continuously while the solutions of equations are converging in newton method.\nyes(y) or no(n)\n");
	fgets(buf3, 256, stdin);
	if(strcmp(buf3, "y\n") == 0 || strcmp(buf3, "yes\n") == 0) {
	  use_lazy_method = 1;
	  break;
	} else if(strcmp(buf3, "n\n") == 0 || strcmp(buf3, "no\n") == 0) {
	  use_lazy_method = 0;
	  break;
	}
	printf("Invalid Input!!!\nSelect yes(y) or no(n)\n");
      }
    }
    if(use_lazy_method == 1) {
      printf("  simulate with lazy mode\n");
    }
    rtn = simulate_implicit(m, &result, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval, &time, order, use_lazy_method, print_amount, mem);
  }

  //write CSV
  if (rtn == NULL) {
    printf("Returned result is NULL\n");
  } else {
    write_csv(rtn, "out.csv"); // for SBML test suite
    /* for more generic simulator
    write_separate_result(rtn,
        "./simulation_results/species_result.dat",
        "./simulation_results/parameter_result.dat",
        "./simulation_results/compartment_result.dat");
     */
  }

  //print result list
  //print_result_list(m, mySp, myParam, myComp);

  //free
  printf("  free all allocated memory\n");
  free_mySBML_objects(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, mem, cp_AST);
  SBMLDocument_free(d);
  return 0;
}
