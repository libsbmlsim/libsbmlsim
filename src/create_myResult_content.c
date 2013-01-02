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

/* create contents of myResult object */
myResult *create_myResult(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], double sim_time, double dt, int print_interval) {
  int i;
  myResult *result;
  int num_of_species = Model_getNumSpecies(m);
  int num_of_parameters = Model_getNumParameters(m);
  int num_of_compartments = Model_getNumCompartments(m);
  int end_cycle = get_end_cycle(sim_time, dt);

  result = (myResult *)malloc(sizeof(myResult));
  result->error_code = NoError;
  result->error_message = NULL;
  result->num_of_columns_sp = num_of_species;
  result->num_of_columns_param = num_of_parameters;
  result->num_of_columns_comp = num_of_compartments;
  result->num_of_rows = end_cycle / print_interval + 1;
  result->column_name_time  = dupstr("time");
  result->column_name_sp    = (const char **)malloc(sizeof(char *) * num_of_species);
  result->column_name_param = (const char **)malloc(sizeof(char *) * num_of_parameters);
  result->column_name_comp  = (const char **)malloc(sizeof(char *) * num_of_compartments);
  result->values_time = (double *)malloc(sizeof(double) * result->num_of_rows);
  result->values_sp = (double *)malloc(sizeof(double) * num_of_species * result->num_of_rows);
  result->values_param = (double *)malloc(sizeof(double) * num_of_parameters * result->num_of_rows);
  result->values_comp = (double *)malloc(sizeof(double) * num_of_compartments * result->num_of_rows);
  for(i=0; i<num_of_species; i++){
    result->column_name_sp[i] = dupstr(Species_getId(mySp[i]->origin));
  }
  for(i=0; i<num_of_parameters; i++){
    result->column_name_param[i] = dupstr(Parameter_getId(myParam[i]->origin));
  }
  for(i=0; i<num_of_compartments; i++){
    result->column_name_comp[i] = dupstr(Compartment_getId(myComp[i]->origin));
  }
  return result;
}

/* create contents of myResult object */
myResult *create_myResultf(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], double sim_time, double dt) {
  int i;
  myResult *result;
  int num_of_species = Model_getNumSpecies(m);
  int num_of_parameters = Model_getNumParameters(m);
  int num_of_compartments = Model_getNumCompartments(m);
  int end_cycle = get_end_cycle(sim_time, dt);

  result = (myResult *)malloc(sizeof(myResult));
  result->error_code = NoError;
  result->error_message = NULL;
  result->num_of_columns_sp = num_of_species;
  result->num_of_columns_param = num_of_parameters;
  result->num_of_columns_comp = num_of_compartments;
  /* new code */
  result->num_of_rows = end_cycle + 1;

  /* new code end */
  result->column_name_time  = dupstr("time");
  result->column_name_sp    = (const char **)malloc(sizeof(char *) * num_of_species);
  result->column_name_param = (const char **)malloc(sizeof(char *) * num_of_parameters);
  result->column_name_comp  = (const char **)malloc(sizeof(char *) * num_of_compartments);
  result->values_time = (double *)malloc(sizeof(double) * result->num_of_rows);
  result->values_sp = (double *)malloc(sizeof(double) * num_of_species * result->num_of_rows);
  result->values_param = (double *)malloc(sizeof(double) * num_of_parameters * result->num_of_rows);
  result->values_comp = (double *)malloc(sizeof(double) * num_of_compartments * result->num_of_rows);
  for(i=0; i<num_of_species; i++){
    result->column_name_sp[i] = dupstr(Species_getId(mySp[i]->origin));
  }
  for(i=0; i<num_of_parameters; i++){
    result->column_name_param[i] = dupstr(Parameter_getId(myParam[i]->origin));
  }
  for(i=0; i<num_of_compartments; i++){
    result->column_name_comp[i] = dupstr(Compartment_getId(myComp[i]->origin));
  }
  return result;
}


/* create myResult object with errorcode and errormessage */
myResult *create_myResult_with_error(LibsbmlsimErrorCode code, const char *message)
{
  myResult *result;
  result = (myResult *)malloc(sizeof(myResult));
  result->error_code = code;
  result->error_message = dupstr(message);

  result->num_of_columns_sp = 0;
  result->num_of_columns_param = 0;
  result->num_of_columns_comp = 0;
  result->num_of_rows = 0;
  result->column_name_time = NULL;
  result->column_name_sp = NULL;
  result->column_name_param = NULL;
  result->column_name_comp = NULL;
  result->values_time = NULL;
  result->values_sp = NULL;
  result->values_param = NULL;
  result->values_comp = NULL;

  return result;
}

/* create myResult object with error and default errormessage */
myResult *create_myResult_with_errorCode(LibsbmlsimErrorCode code)
{
  switch (code) {
    case FileNotFound:
      return create_myResult_with_error(code, "File Not Found");
    case InvalidSBML:
      return create_myResult_with_error(code, "Invalid SBML File");
    case SBMLOperationFailed:
      return create_myResult_with_error(code, "SBML Operation Failed");
    case OutOfMemory:
      return create_myResult_with_error(code, "Out Of Memory");
    case InternalParserError:
      return create_myResult_with_error(code, "Internal SBML Parser Error");
    case SimulationFailed:
      return create_myResult_with_error(code, "Simulation Failed");
    default:
      break;
  }
  return create_myResult_with_error(code, "Unknown Error");
}

SBMLSIM_EXPORT int myResult_isError(myResult *result)
{
  if (result->error_code == NoError)
    return 0;
  return 1;
}

SBMLSIM_EXPORT const char *myResult_getErrorMessage(myResult *result)
{
  return result->error_message;
}

SBMLSIM_EXPORT void __free_myResult(myResult *res)
{
  free_myResult(res);
}

SBMLSIM_EXPORT void free_myResult(myResult *res)
{
  int i;

  if (res->column_name_time != NULL)
    free((void *)res->column_name_time);

  if (res->column_name_sp != NULL) {
    for (i = 0; i < res->num_of_columns_sp; i++)
      free((void *)res->column_name_sp[i]);
    free(res->column_name_sp);
  }

  if (res->column_name_param != NULL) {
    for (i = 0; i < res->num_of_columns_param; i++)
      free((void *)res->column_name_param[i]);
    free(res->column_name_param);
  }
  if (res->column_name_comp != NULL) {
	  for (i = 0; i < res->num_of_columns_comp; i++) {
      free((void *)res->column_name_comp[i]);
	  }
    free(res->column_name_comp);
  }
  if (res->values_time != NULL)
    free(res->values_time);
  if (res->values_sp != NULL)
    free(res->values_sp);
  if (res->values_param != NULL)
    free(res->values_param);
  if (res->values_comp != NULL)
    free(res->values_comp);
  if (res->error_message != NULL)
    free((void *)res->error_message);

  free(res);
}
