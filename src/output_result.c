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
#include <errno.h> 
#include "libsbmlsim/libsbmlsim.h"

SBMLSIM_EXPORT void print_result(myResult* result){
  output_result(result, stdout, ' ');
}

SBMLSIM_EXPORT void write_result(myResult* result, const char* file) {
  char delimiter = ' ';
  print_result_to_file(result, file, delimiter);
}

SBMLSIM_EXPORT void write_csv(myResult* result, const char* file) {
  char delimiter = ',';
  print_result_to_file(result, file, delimiter);
}

FILE* my_fopen(FILE* fp, const char* file, char* mode) {
#ifdef _MSC_VER
  errno_t error;
  char buffer[255];
  if ((error = fopen_s(&fp, file, mode)) != 0) {
    fprintf(stderr, "Failed to open %s: %s\n", file, strerror_s(buffer, sizeof(buffer), error));
    return NULL;
  }
#else
  if ((fp = fopen(file, mode)) == NULL ) {
    fprintf(stderr, "Failed to open %s: %s\n", file, strerror(errno));
    return NULL;
  }
#endif
  return fp;
}

void print_result_to_file(myResult* result, const char* file, char delimiter){
  FILE *fp = NULL;
  if ((fp = my_fopen(fp, file, "w")) != NULL) {
    output_result(result, fp, delimiter);
    fclose(fp);
  }
}

void output_result(myResult* result, FILE* fp, char delimiter){
  int i, j;
  double *value_time_p  = result->values_time;
  double *value_sp_p    = result->values_sp;
  double *value_param_p = result->values_param;
  double *value_comp_p  = result->values_comp;

  /*  Column name (time) */
  fprintf(fp, "%s", result->column_name_time);
  /*  Column name (Species) */
  for (i = 0; i < result->num_of_columns_sp; i++) {
    fprintf(fp, "%c%s", delimiter, result->column_name_sp[i]);
  }
  /*  Column name (Parameters) */
  for (i = 0; i < result->num_of_columns_param; i++) {
    fprintf(fp, "%c%s", delimiter, result->column_name_param[i]);
  }
  /*  Column name (Compartments) */
  for (i = 0; i < result->num_of_columns_comp; i++) {
    fprintf(fp, "%c%s", delimiter, result->column_name_comp[i]);
  }
  fprintf(fp, "\n");

  for (i = 0; i < result->num_of_rows; i++) {
    fprintf(fp, "%.16g", *(value_time_p));
    value_time_p++;
    /*  Species */
    for (j = 0; j < result->num_of_columns_sp; j++) {
      fprintf(fp, "%c%.16g", delimiter, *(value_sp_p));
      value_sp_p++;
    }
    /*  Parameters */
    for (j = 0; j < result->num_of_columns_param; j++) {
      fprintf(fp, "%c%.16g", delimiter, *(value_param_p));
      value_param_p++;
    }
    /*  Compartments */
    for (j = 0; j < result->num_of_columns_comp; j++) {
      fprintf(fp, "%c%.16g", delimiter, *(value_comp_p));
      value_comp_p++;
    }
    fprintf(fp, "\n");
  }
}

SBMLSIM_EXPORT void write_separate_result(myResult* result, const char* file_s, const char* file_p, const char* file_c) {
  FILE *fp_s = NULL;
  FILE *fp_p = NULL;
  FILE *fp_c = NULL;
  int i, j;
  char delimiter = ' ';
  double *value_time_p  = result->values_time;
  double *value_sp_p    = result->values_sp;
  double *value_param_p = result->values_param;
  double *value_comp_p  = result->values_comp;

  if ((fp_s = my_fopen(fp_s, file_s, "w")) == NULL) {
    return;
  }
  if ((fp_p = my_fopen(fp_p, file_p, "w")) == NULL) {
    return;
  }
  if ((fp_c = my_fopen(fp_c, file_c, "w")) == NULL ) {
    return;
  }

  /*  Species */
  for (i = 0; i < result->num_of_rows; i++) {
    fprintf(fp_s, "%.16g", *(value_time_p));
    value_time_p++;
    for (j = 0; j < result->num_of_columns_sp; j++) {
      fprintf(fp_s, "%c%.16g", delimiter, *(value_sp_p));
      value_sp_p++;
    }
    fprintf(fp_s, "\n");
  }
  /*  Parameters */
  value_time_p  = result->values_time;
  for (i = 0; i < result->num_of_rows; i++) {
    fprintf(fp_p, "%.16g", *(value_time_p));
    value_time_p++;
    for (j = 0; j < result->num_of_columns_param; j++) {
      fprintf(fp_p, "%c%.16g", delimiter, *(value_param_p));
      value_param_p++;
    }
    fprintf(fp_p, "\n");
  }
  /*  Compartments */
  value_time_p  = result->values_time;
  for (i = 0; i < result->num_of_rows; i++) {
    fprintf(fp_c, "%.16g", *(value_time_p));
    value_time_p++;
    for (j = 0; j < result->num_of_columns_comp; j++) {
      fprintf(fp_c, "%c%.16g", delimiter, *(value_comp_p));
      value_comp_p++;
    }
    fprintf(fp_c, "\n");
  }
  fclose(fp_s);
  fclose(fp_p);
  fclose(fp_c);
}
