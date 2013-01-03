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
#ifndef LibSBMLSim_MyResult_h
#define LibSBMLSim_MyResult_h

#include "errorcodes.h"

typedef struct myResult{
  LibsbmlsimErrorCode error_code;
  const char *error_message;
  int num_of_rows;
  int num_of_columns_sp;
  int num_of_columns_param;
  int num_of_columns_comp;
  const char *column_name_time;
  const char **column_name_sp;
  const char **column_name_param;
  const char **column_name_comp;
  double *values_time;
  double *values_sp;
  double *values_param;
  double *values_comp;
	/* new code*/
  double* values_time_fordelay;
  int num_of_delay_rows;
} myResult;

#endif /* LibSBMLSim_MyResult_h */
