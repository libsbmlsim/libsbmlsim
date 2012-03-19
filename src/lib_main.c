#include "libsbmlsim/libsbmlsim.h"

/* libSBMLSimulator API */

SBMLSIM_EXPORT myResult* simulateSBMLFromString(const char* str, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method) {
  SBMLDocument_t* d;
  Model_t* m;
  myResult *rtn;
  d = readSBMLFromString(str);
  m = SBMLDocument_getModel(d);
  printf("simulation start\n");
  rtn = simulateSBMLModel(m, sim_time, dt, print_interval, print_amount, method, use_lazy_method);
  printf("simulation end\n");
  SBMLDocument_free(d);
  return rtn;
}

SBMLSIM_EXPORT myResult* simulateSBMLModel(Model_t *m, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method){
  double time = 0;
  int order = 0;
  int is_explicit = 0;
  char *method_name;
  unsigned int num_of_species;
  unsigned int num_of_parameters;
  unsigned int num_of_compartments;
  unsigned int num_of_reactions;
  unsigned int num_of_rules;
  unsigned int num_of_events;
  unsigned int num_of_initialAssignments;
  mySpecies **mySp;
  myParameter **myParam;
  myCompartment **myComp;
  myReaction **myRe;
  myRule **myRu;
  myEvent **myEv;
  myInitialAssignment **myInitAssign;
  /* prepare myAlgebraicEquations */
  myAlgebraicEquations *myAlgEq = NULL;
  /* prepare timeVariantAssignments */
  timeVariantAssignments *timeVarAssign = NULL;
  /* prepare return value */
  myResult *result, *rtn;

  allocated_memory *mem;
  copied_AST *cp_AST;

  mem = (allocated_memory*)malloc(sizeof(allocated_memory));
  mem->num_of_allocated_memory = 0;

  cp_AST = (copied_AST*)malloc(sizeof(copied_AST));
  cp_AST->num_of_copied_AST = 0;

  /* prepare mySpecies */
  num_of_species = Model_getNumSpecies(m);
  /* mySpecies *mySp[num_of_species]; */
  mySp = (mySpecies**)malloc(sizeof(mySpecies*) * num_of_species);
  /* prepare myParameters */
  num_of_parameters = Model_getNumParameters(m);
  /* myParameter *myParam[num_of_parameters]; */
  myParam = (myParameter**)malloc(sizeof(myParameter*) * num_of_parameters);
  /* prepare myCompartments */
  num_of_compartments = Model_getNumCompartments(m);
  /* myCompartment *myComp[num_of_compartments]; */
  myComp = (myCompartment**)malloc(sizeof(mySpecies*) * num_of_compartments);
  /* prepare myReactions */
  num_of_reactions = Model_getNumReactions(m);
  /* myReaction *myRe[num_of_reactions]; */
  myRe = (myReaction**)malloc(sizeof(myReaction*) * num_of_reactions);
  /* prepare myRules */
  num_of_rules = Model_getNumRules(m);
  /* myRule *myRu[num_of_rules]; */
  myRu = (myRule**)malloc(sizeof(myRule*) * num_of_rules);
  /* prepare myEvents */
  num_of_events = Model_getNumEvents(m);  
  /* myEvent *myEv[num_of_events]; */
  myEv = (myEvent**)malloc(sizeof(myEvent*) * num_of_events);
  /* prepare myInitial Assignments */
  num_of_initialAssignments = Model_getNumInitialAssignments(m);
  /* myInitialAssignment *myInitAssign[num_of_initialAssignments]; */
  myInitAssign = (myInitialAssignment**)malloc(sizeof(myInitialAssignment*) * num_of_initialAssignments);
  /* create myObjects */
  create_mySBML_objects(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, &myAlgEq, &timeVarAssign, sim_time, dt, &time, mem, cp_AST);
  /* create myResult */
  result = create_myResult(m, mySp, myParam, myComp, sim_time, dt, print_interval);

  switch(method) {
    case MTHD_RUNGE_KUTTA: /*  Runge-Kutta */
      method_name = MTHD_NAME_RUNGE_KUTTA;
      break;
    case MTHD_BACKWARD_EULER: /*  Backward-Euler */
      method_name = MTHD_NAME_BACKWARD_EULER;
      break;
    case MTHD_CRANK_NICOLSON: /*  Crank-Nicolson (Adams-Moulton 2) */
      method_name = MTHD_NAME_CRANK_NICOLSON;
      break;
    case MTHD_ADAMS_MOULTON_3: /*  Adams-Moulton 3 */
      method_name = MTHD_NAME_ADAMS_MOULTON_3;
      break;
    case MTHD_ADAMS_MOULTON_4: /*  Adams-Moulton 4 */
      method_name = MTHD_NAME_ADAMS_MOULTON_4;
      break;
    case MTHD_BACKWARD_DIFFERENTIATION_2: /*  Backward-Differentiation 2 */
      method_name = MTHD_NAME_BACKWARD_DIFFERENTIATION_2;
      break;
    case MTHD_BACKWARD_DIFFERENTIATION_3: /*  Backward-Differentiation 3 */
      method_name = MTHD_NAME_BACKWARD_DIFFERENTIATION_3;
      break;
    case MTHD_BACKWARD_DIFFERENTIATION_4: /*  Backward-Differentiation 4 */
      method_name = MTHD_NAME_BACKWARD_DIFFERENTIATION_4;
      break;
    case MTHD_EULER: /*  Euler (Adams-Bashforth) */
      method_name = MTHD_NAME_EULER;
      break;
    case MTHD_ADAMS_BASHFORTH_2: /*  Adams-Bashforth 2 */
      method_name = MTHD_NAME_ADAMS_BASHFORTH_2;
      break;
    case MTHD_ADAMS_BASHFORTH_3: /*  Adams-Bashforth 3 */
      method_name = MTHD_NAME_ADAMS_BASHFORTH_3;
      break;
    case MTHD_ADAMS_BASHFORTH_4: /*  Adams-Bashforth 4 */
      method_name = MTHD_NAME_ADAMS_BASHFORTH_4;
      break;
    default:
      method = MTHD_RUNGE_KUTTA;
      method_name = MTHD_NAME_RUNGE_KUTTA;
      break;
  }
  order = method / 10;
  is_explicit = method % 10;
  TRACE(("simulate with %s\n", method_name));

  /* simulation */
  if (is_explicit == 1) {
    rtn = simulate_explicit(m, result, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval, &time, order, print_amount, mem);
  }else{
    rtn = simulate_implicit(m, result, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval, &time, order, use_lazy_method, print_amount, mem);
  }

  /* free */
  free_mySBML_objects(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, mem, cp_AST);
  if (rtn == NULL)
    free_myResult(result);
  return rtn;
}
