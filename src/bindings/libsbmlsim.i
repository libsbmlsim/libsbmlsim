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
/* libsbmlsim.i */
%module libsbmlsim

%{
#include "../../src/libsbmlsim/myResult.h"
extern myResult* simulateSBMLFromFile(const char *file, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method);
extern myResult* simulateSBMLFromString(const char *str, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method);
extern void print_result(myResult* result);
extern void write_result(myResult* result, char* file);
extern void write_csv(myResult* result, char* file);
extern void write_separate_result(myResult* result, char* file_s, char* file_p, char* file_c);
extern void __free_myResult(myResult *result);
typedef int BOOLEAN;
%}

%include "../../src/libsbmlsim/methods.h"
%include "../../src/libsbmlsim/errorcodes.h"

#ifdef SWIGJAVA
%apply bool { BOOLEAN };
#endif

#ifdef SWIGPYTHON
%typemap(out) BOOLEAN %{
  if ($1)
    $result = Py_True;
  else
    $result = Py_False;
%}
#endif

#ifdef SWIGRUBY
%typemap(out) BOOLEAN %{
  if ($1)
    $result = Qtrue;
  else
    $result = Qfalse;
%}
#endif

#ifdef SWIGPERL
%typemap(out) BOOLEAN %{
  $result = boolSV($1);
%}
#endif

/* %include "src/libsbmlsim/myResult.h" */
typedef struct myResult {
  LibsbmlsimErrorCode error_code;
%immutable;
  char *error_message;
%mutable;
  int num_of_rows;
  int num_of_columns_sp;
  int num_of_columns_param;
  int num_of_columns_comp;
%immutable;
  char *column_name_time;
  char **column_name_sp;
  char **column_name_param;
  char **column_name_comp;
%mutable;
  double *values_time;
  double *values_sp;
  double *values_param;
  double *values_comp;
  double *values_time_fordelay;
  int num_of_delay_rows;
} myResult;

%newobject simulateSBMLFromFile;
%newobject simulateSBMLFromString;
extern myResult* simulateSBMLFromFile(const char *file, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method);
extern myResult* simulateSBMLFromString(const char *str, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method);
extern void print_result(myResult* result);
extern void write_result(myResult* result, char* file);
extern void write_csv(myResult* result, char* file);
extern void write_separate_result(myResult* result, char* file_s, char* file_p, char* file_c);

%extend myResult {
  myResult() {
    myResult *res;
    res = (myResult *)malloc(sizeof(myResult));
    return res;
  }

  ~myResult() {
    __free_myResult($self);
  }

  BOOLEAN isError() {
    return $self->error_code != NoError;
  }

  LibsbmlsimErrorCode getErrorCode() {
    return $self->error_code;
  }

  const char *getErrorMessage() {
    return $self->error_message;
  }

  int getNumOfRows() {
    return $self->num_of_rows;
  }

  int getNumOfSpecies() {
    return $self->num_of_columns_sp;
  }

  int getNumOfParameters() {
    return $self->num_of_columns_param;
  }

  int getNumOfCompartments() {
    return $self->num_of_columns_comp;
  }

  const char *getTimeName() {
    return $self->column_name_time;
  }

  const char *getSpeciesNameAtIndex(int index) {
    if (index < 0 || index >= $self->num_of_columns_sp || $self->error_code != NoError)
      return NULL;
    return $self->column_name_sp[index];
  }

  const char *getParameterNameAtIndex(int index) {
    if (index < 0 || index >= $self->num_of_columns_param || $self->error_code != NoError)
      return NULL;
    return $self->column_name_param[index];
  }

  const char *getCompartmentNameAtIndex(int index) {
    if (index < 0 || index >= $self->num_of_columns_comp || $self->error_code != NoError)
      return NULL;
    return $self->column_name_comp[index];
  }

  double getTimeValueAtIndex(int index) {
    if (index < 0 || index >= $self->num_of_rows || $self->error_code != NoError)
      return -0.0;
    return $self->values_time[index];
  }

  double getSpeciesValueAtIndex(char *sname, int index) {
    int i, spindex;
    spindex = -1;
    if (index < 0 || index >= $self->num_of_rows || $self->error_code != NoError)
      return -0.0;
    for (i = 0; i < $self->num_of_columns_sp; i++) {
      if (strcmp($self->column_name_sp[i], sname) == 0) {
        spindex = i;
        break;
      }
    }
    if (spindex == -1)
      return -0.0;
    return $self->values_sp[index * $self->num_of_columns_sp + spindex];
  }

  double getParameterValueAtIndex(char *pname, int index) {
    int i, pindex;
    pindex = -1;
    if (index < 0 || index >= $self->num_of_rows || $self->error_code != NoError)
      return -0.0;
    for (i = 0; i < $self->num_of_columns_param; i++) {
      if (strcmp($self->column_name_param[i], pname) == 0) {
        pindex = i;
        break;
      }
    }
    if (pindex == -1)
      return -0.0;
    return $self->values_param[index];
  }

  double getCompartmentValueAtIndex(char *cname, int index) {
    int i, cindex;
    cindex = -1;
    if (index < 0 || index >= $self->num_of_rows || $self->error_code != NoError)
      return -0.0;
    for (i = 0; i < $self->num_of_columns_comp; i++) {
      if (strcmp($self->column_name_comp[i], cname) == 0) {
        cindex = i;
        break;
      }
    }
    if (cindex == -1)
      return -0.0;
    return $self->values_comp[index];
  }

};
