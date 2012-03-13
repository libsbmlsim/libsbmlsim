#include "libsbmlsim/libsbmlsim.h"

/* libSBMLSimulator API */

myResult* simulateSBMLFromString(const char* str, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method) {
  SBMLDocument_t* d;
  Model_t* m;
  myResult result;
  myResult *rtn;
  d = readSBMLFromString(str);
  m = SBMLDocument_getModel(d);
  rtn = simulateSBMLModel(m, &result, sim_time, dt, print_interval, print_amount, method, use_lazy_method);
  SBMLDocument_free(d);
  return rtn;
}

myResult* simulateSBMLModel(Model_t *m, myResult* result, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method){
  int is_explicit = 0;
  double time = 0;
  int order = 0;

  allocated_memory *mem;
  mem = (allocated_memory*)malloc(sizeof(allocated_memory));
  mem->num_of_allocated_memory = 0;

  copied_AST *cp_AST;
  cp_AST = (copied_AST*)malloc(sizeof(copied_AST));
  cp_AST->num_of_copied_AST = 0;

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
  //create myResult
  create_myResult_content(m, result, mySp, myParam, myComp, sim_time, dt, print_interval);
  //prepare return value
  myResult* rtn;

  switch(method) {
    case 1:
      dbg_printf("simulate with runge kutta\n");
      order = 4;
      is_explicit = 1;
      break;
    case 2:
      dbg_printf("simulate with Backward-Euler\n");
      order = 0;
      break;
    case 3:
      dbg_printf("simulate with AM2(Crank-Nicolson)\n");
      order = 1;
      break;
    case 4:
      dbg_printf("simulate with AM3\n");
      order = 2;
      break;
    case 5:
      dbg_printf("simulate with AM4\n");
      order = 3;
      break;
    case 6:
      dbg_printf("simulate with BD2\n");
      order = 4;
      break;
    case 7:
      dbg_printf("simulate with BD3\n");
      order = 5;
      break;
    case 8:
      dbg_printf("simulate with BD4\n");
      order = 6;
      break;
    case 9:
      dbg_printf("simulate with AB1(Euler)\n");
      order = 0;
      is_explicit = 1;
      break;
    case 10:
      dbg_printf("simulate with AB2\n");
      order = 1;
      is_explicit = 1;
      break;
    case 11:
      dbg_printf("simulate with AB3\n");
      order = 2;
      is_explicit = 1;
      break;
    case 12:
      dbg_printf("simulate with AB4\n");
      order = 3;
      is_explicit = 1;
      break;
    default:
      dbg_printf("simulate with runge kutta\n");
      order = 4;
      break;
  }

  //simulation
  if (is_explicit == 1) {
    rtn = simulate_explicit(m, result, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval, &time, order, print_amount, mem);
  }else{
    rtn = simulate_implicit(m, result, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval, &time, order, use_lazy_method, print_amount, mem);
  }

  //print result list
  //print_result_list(m, mySp, myParam, myComp);

  //free
  free_mySBML_objects(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, mem, cp_AST);
  return rtn;
}
