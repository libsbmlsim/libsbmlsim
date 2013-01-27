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

/* libSBMLSimulator API */

SBMLSIM_EXPORT myResult* simulateSBMLFromFile(const char* file, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method) {
  SBMLDocument_t* d;
  Model_t* m;
  myResult *rtn;
  unsigned int err_num;
  double atol = 0.0;
  double rtol = 0.0;
  double facmax = 0.0;
  d = readSBMLFromFile(file);
  if (d == NULL)
    return create_myResult_with_errorCode(Unknown);
  err_num = SBMLDocument_getNumErrors(d);
  if (err_num > 0) {
    const XMLError_t *err = (const XMLError_t *)SBMLDocument_getError(d, 0);
    if (XMLError_isError(err) || XMLError_isFatal(err)) {
      XMLErrorCode_t errcode = XMLError_getErrorId(err);
      switch (errcode) {
        case XMLFileUnreadable:
          rtn = create_myResult_with_errorCode(FileNotFound);
          break;
        case XMLFileUnwritable:
        case XMLFileOperationError:
        case XMLNetworkAccessError:
          rtn = create_myResult_with_errorCode(SBMLOperationFailed);
          break;
        case InternalXMLParserError:
        case UnrecognizedXMLParserCode:
        case XMLTranscoderError:
          rtn = create_myResult_with_errorCode(InternalParserError);
          break;
        case XMLOutOfMemory:
          rtn = create_myResult_with_errorCode(OutOfMemory);
          break;
        case XMLUnknownError:
          rtn = create_myResult_with_errorCode(Unknown);
          break;
        default:
          rtn = create_myResult_with_errorCode(InvalidSBML);
          break;
      }
      SBMLDocument_free(d);
      return rtn;
    }
  }
  m = SBMLDocument_getModel(d);
  rtn = simulateSBMLModel(m, sim_time, dt, print_interval, print_amount, method, use_lazy_method, atol, rtol, facmax);
  if (rtn == NULL)
    rtn = create_myResult_with_errorCode(SimulationFailed);
  SBMLDocument_free(d);
  return rtn;
}

SBMLSIM_EXPORT myResult* simulateSBMLFromString(const char* str, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method) {
  SBMLDocument_t* d;
  Model_t* m;
  myResult *rtn;
  unsigned int err_num;
  double atol = 0.0;
  double rtol = 0.0;
  double facmax = 0.0;
  d = readSBMLFromString(str);
  if (d == NULL)
    return create_myResult_with_errorCode(Unknown);
  err_num = SBMLDocument_getNumErrors(d);
  if (err_num > 0) {
    const XMLError_t *err = (const XMLError_t *)SBMLDocument_getError(d, 0);
    if (XMLError_isError(err) || XMLError_isFatal(err)) {
      XMLErrorCode_t errcode = XMLError_getErrorId(err);
      switch (errcode) {
        case XMLFileUnreadable:
          rtn = create_myResult_with_errorCode(FileNotFound);
          break;
        case XMLFileUnwritable:
        case XMLFileOperationError:
        case XMLNetworkAccessError:
          rtn = create_myResult_with_errorCode(SBMLOperationFailed);
          break;
        case InternalXMLParserError:
        case UnrecognizedXMLParserCode:
        case XMLTranscoderError:
          rtn = create_myResult_with_errorCode(InternalParserError);
          break;
        case XMLOutOfMemory:
          rtn = create_myResult_with_errorCode(OutOfMemory);
          break;
        case XMLUnknownError:
          rtn = create_myResult_with_errorCode(Unknown);
          break;
        default:
          rtn = create_myResult_with_errorCode(InvalidSBML);
          break;
      }
      SBMLDocument_free(d);
      return rtn;
    }
  }
  m = SBMLDocument_getModel(d);
  rtn = simulateSBMLModel(m, sim_time, dt, print_interval, print_amount, method, use_lazy_method, atol, rtol, facmax);
  if (rtn == NULL)
    rtn = create_myResult_with_errorCode(SimulationFailed);
  SBMLDocument_free(d);
  return rtn;
}

SBMLSIM_EXPORT myResult* simulateSBMLModel(Model_t *m, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method, double atol, double rtol, double facmax){
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
  myResult *result, *rtn = NULL;
  /* Variables for bifurcation analysis */
  char buf1[256], buf2[256], buf3[256], buf4[256], buf5[256], buf6[256];
  boolean use_bifurcation_analysis = false; /* Bifurcation is turned off for this release */
  boolean state_variable_exists = false;
  boolean bifurcation_parameter_exists = false;
  boolean bif_param_is_local = false;
  boolean is_variable_step = false;
  char* sta_var_id = NULL;
  char* bif_param_id = NULL;

  double bif_param_min = 0;
  double bif_param_max = 0;
  double bif_param_stepsize = 0.01;

  /* to exclude the transition state */
  double transition_time = 0;

  char* tmp;
  unsigned int i, a, b;
  
  /* for variable stepsize */
  int err_zero_flag = 0;

  allocated_memory *mem;
  copied_AST *cp_AST;

  mem = allocated_memory_create();
  cp_AST = copied_AST_create();

  /* Check atol, rtol and facmax, whether it is set to 0.0 */
  if (atol == 0.0) {
    atol = ABSOLUTE_ERROR_TOLERANCE;
  }
  if (rtol == 0.0) {
    rtol = RELATIVE_ERROR_TOLERANCE;
  }
  if (facmax == 0.0) {
    facmax = DEFAULT_FACMAX;
  }
  
  /* determin bifurcation analysis condition */
  if (use_bifurcation_analysis) {
	  while(1) {
		  printf("--- Bifurcation Analysis mode ---\n");
		  printf("[1] select 1 state variable(species ID) and 1 bifurcation parameter(ID) from the following list\n");
		  /*show the list of species and parameter*/
		  show_sp(m);
		  printf("Species ID : ");
		  tmp = fgets(buf1, 256, stdin);
		  chomp(buf1);
		  for(i = 0; i < Model_getNumSpecies(m); i++) {
			  if (strcmp(buf1, Species_getId((Species_t*)ListOf_get(Model_getListOfSpecies(m), i))) == 0) {
				  printf("------ %s is selected. \n", buf1);
				  state_variable_exists = true;
				  break;
			  }
		  }
		  show_para(m);
		  printf("Parameter ID : ");
		  tmp = fgets(buf2, 256, stdin);
		  chomp(buf2);
		  for(i = 0; i < Model_getNumParameters(m); i++) {
			  if (strcmp(buf2, Parameter_getId((Parameter_t*)ListOf_get(Model_getListOfParameters(m), i))) == 0) {
				  printf("------ %s is selected. \n", buf2);
				  bif_param_is_local = false;
				  bifurcation_parameter_exists = true;
				  break;
			  }
		  }
		  if (bifurcation_parameter_exists == false) {
			  for(a = 0; a < Model_getNumReactions(m); a++) {
				  for(b = 0; b < ListOf_size(KineticLaw_getListOfParameters(Reaction_getKineticLaw((Reaction_t*)ListOf_get(Model_getListOfReactions(m), a)))); b++){
					  if (strcmp(buf2, Parameter_getId(KineticLaw_getParameter(Reaction_getKineticLaw((Reaction_t*)ListOf_get(Model_getListOfReactions(m), a) ), b))) == 0) {
						  printf("------ %s is selected. \n", buf2);
						  bifurcation_parameter_exists = true;
						  bif_param_is_local = true;
						  break;
					  }
				  }
			  }
		  }
		  if(!state_variable_exists) {
			  printf("[ERROR] Species : %s  is not found.\n", buf1);
		  }
		  if(!bifurcation_parameter_exists) {
			  printf("[ERROR] Paraneter : %s is not found.\n", buf2);
		  }
		  if(state_variable_exists == true && bifurcation_parameter_exists == true) {
			  while(1) {
				  printf("[2] set the transition duration (less than simulation time)\n");
				  tmp = fgets(buf3, 256, stdin);
				  chomp(buf3);
				  if(str_is_number(buf3) == 1 && atof(buf3) < sim_time) {
					  break;
				  }
				  printf("is not a number or not less than simulation time! \n");
			  }
			  while(1) {
				  printf("[3] specify the value range of the parameter you selected.\n");
				  while(1) {
					  printf("the minimum value : ");
					  tmp = fgets(buf4, 256, stdin);
					  chomp(buf4);
					  if(str_is_number(buf4)){
						  break;
					  }
					  printf("  is not a number!\n");
				  }
				  while(1) {
					  printf("the maximum value : ");
					  tmp = fgets(buf5, 256, stdin);
					  chomp(buf5);
					  if(str_is_number(buf5)){
						  break;
					  }
					  printf("  is not a number!\n");
				  }
				  if (atof(buf4) < atof(buf5)) {
					  break;
				  } else {
					  printf("\nset two values so that the maximum value is larger than the minimum one\n");
				  }
			  }
			  while(1) {
				  printf("the step size of bifurcation parameter : ");
				  tmp = fgets(buf6, 256, stdin);
				  chomp(buf6);
				  if(str_is_number(buf6)){
					  break;
					  }
				  printf("  is not a number!\n");
			  }
			  sta_var_id = (char*)malloc(sizeof(char) * strlen(buf1));
			  strcpy(sta_var_id, buf1);
			  bif_param_id = (char*)malloc(sizeof(char) * strlen(buf2));
			  strcpy(bif_param_id, buf2);
			  sscanf(buf3, "%lf", &transition_time);
			  sscanf(buf4, "%lf", &bif_param_min);
			  sscanf(buf5, "%lf", &bif_param_max);
			  sscanf(buf6, "%lf", &bif_param_stepsize);
			  break;
		  }
	  }
  }

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
    case MTHD_BACKWARD_DIFFERENCE_2: /*  Backward-Difference 2 */
      method_name = MTHD_NAME_BACKWARD_DIFFERENCE_2;
      break;
    case MTHD_BACKWARD_DIFFERENCE_3: /*  Backward-Difference 3 */
      method_name = MTHD_NAME_BACKWARD_DIFFERENCE_3;
      break;
    case MTHD_BACKWARD_DIFFERENCE_4: /*  Backward-Difference 4 */
      method_name = MTHD_NAME_BACKWARD_DIFFERENCE_4;
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
    /* Variable Step Size */
    case MTHD_RUNGE_KUTTA_FEHLBERG_5: /*  Runge-Kutta-Fehlberg */
      method_name = MTHD_NAME_RUNGE_KUTTA_FEHLBERG_5;
      is_variable_step = true;
      break;
    case MTHD_CASH_KARP: /*  Cash-Karp */
      method_name = MTHD_NAME_CASH_KARP;
      is_variable_step = true;
      break;
    default:
      method = MTHD_RUNGE_KUTTA;
      method_name = MTHD_NAME_RUNGE_KUTTA;
      break;
  }
  order = method / 10;
  is_explicit = method % 10;
  TRACE(("simulate with %s\n", method_name));

  /* create myObjects */
  if (is_variable_step) {
    create_mySBML_objectsf(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, &myAlgEq, &timeVarAssign, sim_time, dt, &time, mem, cp_AST, print_interval);
  } else {
    create_mySBML_objects(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, &myAlgEq, &timeVarAssign, sim_time, dt, &time, mem, cp_AST);
  }

  /* create myResult */
  if (is_variable_step) {
    result = create_myResultf(m, mySp, myParam, myComp, sim_time, dt);
  } else {
    result = create_myResult(m, mySp, myParam, myComp, sim_time, dt, print_interval);
  }

  /* simulation */
  if (is_variable_step) {
    /* if (order == 5 || order == 6) { */
    if (is_explicit) {
      rtn = simulate_explicitf(m, result, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval, &time, order, print_amount, mem, atol, rtol, facmax, cp_AST, &err_zero_flag);
    }
  } else {  /* Fixed step size */
    if (is_explicit) {
    rtn = simulate_explicit(m, result, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval, &time, order, print_amount, mem);
  }else{
    rtn = simulate_implicit(m, result, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval, &time, order, use_lazy_method, print_amount, mem);
  }
  }

  /* bifurcation analysis */
  if(use_bifurcation_analysis) {
	  rtn = bifurcation_analysis(m, sim_time, dt, print_interval, time, order, print_amount, use_lazy_method, is_explicit, num_of_species, num_of_parameters, num_of_compartments, num_of_reactions, num_of_rules, num_of_events, num_of_initialAssignments, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, mem, cp_AST, result, rtn, bif_param_is_local, sta_var_id, bif_param_id, bif_param_min, bif_param_max, bif_param_stepsize, transition_time);
  }
  /* after bifurcation analysis */
  if (use_bifurcation_analysis) {
	  free(sta_var_id);
	  free(bif_param_id);
  }
  /* free */
  if (use_bifurcation_analysis == 0 || (use_bifurcation_analysis == 1 && bif_param_is_local == false)) {
	  free_mySBML_objects(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, mem, cp_AST);
  }

  if (rtn == NULL)
    free_myResult(result);
  return rtn;
}

SBMLSIM_EXPORT myResult* simulateSBMLModelf(Model_t *m, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method, double atol, double rtol, double facmax, int use_bifurcation_analysis){
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
  myResult *result, *rtn = NULL;
  /* Variables for bifurcation analysis */
  char buf1[256], buf2[256], buf3[256], buf4[256], buf5[256], buf6[256];
  boolean state_variable_exists = false;
  boolean bifurcation_parameter_exists = false;
  boolean bif_param_is_local = false;
  char* sta_var_id = NULL;
  char* bif_param_id = NULL;
  double bif_param_min = 0;
  double bif_param_max = 0;
  double bif_param_stepsize = 0.01;
  char* tmp;
  unsigned int i, a, b;
  /* to exclude the transition state */
  double transition_time = 0;
  /* for variable stepsize */
  int err_zero_flag = 0;
  

  allocated_memory *mem;
  copied_AST *cp_AST;
  
  /* This code is not called from library, but just avoid warnigns,
   * we will asign a value to unused variable
   */
  use_lazy_method = 0;

  mem = allocated_memory_create();
  cp_AST = copied_AST_create();

  /* determin bifurcation analysis condition*/
  if (use_bifurcation_analysis) {
	  while(1) {
		  printf("--- Bifurcation Analysis mode ---\n");
		  printf("[1] select 1 state variable(species ID) and 1 bifurcation parameter(ID) from the following list\n");
		  /*show the list of species and parameter*/
		  show_sp(m);
		  printf("Species ID : ");
		  tmp = fgets(buf1, 256, stdin);
		  chomp(buf1);
		  for(i = 0; i < Model_getNumSpecies(m); i++) {
			  if (strcmp(buf1, Species_getId((Species_t*)ListOf_get(Model_getListOfSpecies(m), i))) == 0) {
				  printf("------ %s is selected. \n", buf1);
				  state_variable_exists = true;
				  break;
			  }
		  }
		  show_para(m);
		  printf("Parameter ID : ");
		  tmp = fgets(buf2, 256, stdin);
		  chomp(buf2);
		  for(i = 0; i < Model_getNumParameters(m); i++) {
			  if (strcmp(buf2, Parameter_getId((Parameter_t*)ListOf_get(Model_getListOfParameters(m), i))) == 0) {
				  printf("------ %s is selected. \n", buf2);
				  bif_param_is_local = false;
				  bifurcation_parameter_exists = true;
				  break;
			  }
		  }
		  if (bifurcation_parameter_exists == false) {
			  for(a = 0; a < Model_getNumReactions(m); a++) {
				  for(b = 0; b < ListOf_size(KineticLaw_getListOfParameters(Reaction_getKineticLaw((Reaction_t*)ListOf_get(Model_getListOfReactions(m), a)))); b++){
					  if (strcmp(buf2, Parameter_getId(KineticLaw_getParameter(Reaction_getKineticLaw((Reaction_t*)ListOf_get(Model_getListOfReactions(m), a) ), b))) == 0) {
						  printf("------ %s is selected. \n", buf2);
						  bifurcation_parameter_exists = true;
						  bif_param_is_local = true;
						  break;
					  }
				  }
			  }
		  }
		  if(!state_variable_exists) {
			  printf("[ERROR] Species : %s  is not found.\n", buf1);
		  }
		  if(!bifurcation_parameter_exists) {
			  printf("[ERROR] Paraneter : %s is not found.\n", buf2);
		  }
		  if(state_variable_exists == true && bifurcation_parameter_exists == true) {
			  while(1) {
				  printf("[2] set the transition duration (less than simulation time)\n");
				  tmp = fgets(buf3, 256, stdin);
				  chomp(buf3);
				  if(str_is_number(buf3) == 1 && atof(buf3) < sim_time) {
					  break;
				  }
				  printf("is not a number or not less than simulation time! \n");
			  }
			  while(1) {
				  printf("[3] specify the value range of the parameter you selected.\n");
				  while(1) {
					  printf("the minimum value : ");
					  tmp = fgets(buf4, 256, stdin);
					  chomp(buf4);
					  if(str_is_number(buf4)){
						  break;
					  }
					  printf("  is not a number!\n");
				  }
				  while(1) {
					  printf("the maximum value : ");
					  tmp = fgets(buf5, 256, stdin);
					  chomp(buf5);
					  if(str_is_number(buf5)){
						  break;
					  }
					  printf("  is not a number!\n");
				  }
				  if (atof(buf4) < atof(buf5)) {
					  break;
				  } else {
					  printf("\nset two values so that the maximum value is larger than the minimum one\n");
				  }
			  }
			  while(1) {
				  printf("the step size of bifurcation parameter : ");
				  tmp = fgets(buf6, 256, stdin);
				  chomp(buf6);
				  if(str_is_number(buf6)){
					  break;
					  }
				  printf("  is not a number!\n");
			  }
			  sta_var_id = (char*)malloc(sizeof(char) * strlen(buf1));
			  strcpy(sta_var_id, buf1);
			  bif_param_id = (char*)malloc(sizeof(char) * strlen(buf2));
			  strcpy(bif_param_id, buf2);
			  sscanf(buf3, "%lf", &transition_time);
			  sscanf(buf4, "%lf", &bif_param_min);
			  sscanf(buf5, "%lf", &bif_param_max);
			  sscanf(buf6, "%lf", &bif_param_stepsize);
			  break;
		  }
	  }
  }

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
  create_mySBML_objectsf(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, &myAlgEq, &timeVarAssign, sim_time, dt, &time, mem, cp_AST, print_interval);
  /* create myResult */
  result = create_myResultf(m, mySp, myParam, myComp, sim_time, dt);

  switch(method) {
    case MTHD_RUNGE_KUTTA_FEHLBERG_5: /*  Runge-Kutta-Fehlberg */
      method_name = MTHD_NAME_RUNGE_KUTTA_FEHLBERG_5;
      break;
    case MTHD_CASH_KARP: /*  Cash-Karp */
      method_name = MTHD_NAME_CASH_KARP;
      break;
    default:
      fprintf(stderr, "implicit numerical method\n");
      exit(1);
      break;
  }
  order = method / 10;
  is_explicit = method % 10;
  TRACE(("simulate with %s\n", method_name));
  /* simulation */
  if (is_explicit == 1) {
	  /* Runge-Kutta Fehlberg, variable step-size numerical integration */
	  if (order == 5 || order == 6) {
 		  rtn = simulate_explicitf(m, result, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval, &time, order, print_amount, mem, atol, rtol, facmax, cp_AST, &err_zero_flag);
		  //write_csv(rtn, "variable.csv"); /* save the result data by variable step-size */
		  if (err_zero_flag == 0) {
			  //rtn = approximate_result_linearly(m, rtn, mySp, myParam, myComp, sim_time, dt);
		  }
	  }
  }

  /*bifurcation analysis (XXX must be fixed)*/
/*   if(use_bifurcation_analysis) { */
/* 	  rtn = bifurcation_analysis(m, sim_time, dt, print_interval, time, order, print_amount, use_lazy_method, is_explicit, num_of_species, num_of_parameters, num_of_compartments, num_of_reactions, num_of_rules, num_of_events, num_of_initialAssignments, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, mem, cp_AST, result, rtn, bif_param_is_local, sta_var_id, bif_param_id, bif_param_min, bif_param_max, bif_param_stepsize, transition_time); */
/*   } */

  /* after bifurcation analysis*/
  if (use_bifurcation_analysis) {
	  free(sta_var_id);
	  free(bif_param_id);
  }
  /* free */
  if (use_bifurcation_analysis == 0 || (use_bifurcation_analysis == 1 && bif_param_is_local == false)) {
  free_mySBML_objects(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, mem, cp_AST);
  }
  if (rtn == NULL)
    free_myResult(result);
  return rtn;
}
