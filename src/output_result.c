#include <errno.h> 
#include "libsbmlsim/libsbmlsim.h"

void print_result(myResult* result){
  output_result(result, stdout, ' ');
}

void write_result(myResult* result, char* file) {
  char delimiter = ' ';
  print_result_to_file(result, file, delimiter);
}

void write_csv(myResult* result, char* file) {
  char delimiter = ',';
  print_result_to_file(result, file, delimiter);
}

void print_result_to_file(myResult* result, char* file, char delimiter){
  FILE *fp;

  fp = fopen(file, "w");
  if (fp == NULL) {
    fprintf(stderr, "File %s open failed.\n", file);
  } else {
    output_result(result, fp, delimiter);
    fclose(fp);
  }
}

void output_result0(myResult* result, FILE* fp, char delimiter){
  int i, j;
  double *value_p = result->values;

  for (i = 0; i < result->num_of_columns; i++) {
    if (i == 0) {
      fprintf(fp, "%s", result->column_name[i]);
    } else {
      fprintf(fp, "%c%s", delimiter, result->column_name[i]);
    }
  }
  fprintf(fp, "\n");
  for (i = 0; i < result->num_of_rows; i++) {
    for (j = 0; j < result->num_of_columns; j++) {
      if (j == 0) {
        fprintf(fp, "%.16g", *(value_p));
      } else {
        fprintf(fp, "%c%.16g", delimiter, *(value_p));
      }
      value_p++;
    }
    fprintf(fp, "\n");
  }
}

void output_result(myResult* result, FILE* fp, char delimiter){
  int i, j;
  double *value_time_p  = result->values_time;
  double *value_sp_p    = result->values_sp;
  double *value_param_p = result->values_param;
  double *value_comp_p  = result->values_comp;

  // Column name (time)
  fprintf(fp, "%s", result->column_name_time);
  // Column name (Species)
  for (i = 0; i < result->num_of_columns_sp; i++) {
    fprintf(fp, "%c%s", delimiter, result->column_name_sp[i]);
  }
  // Column name (Parameters)
  for (i = 0; i < result->num_of_columns_param; i++) {
    fprintf(fp, "%c%s", delimiter, result->column_name_param[i]);
  }
  // Column name (Compartments)
  for (i = 0; i < result->num_of_columns_comp; i++) {
    fprintf(fp, "%c%s", delimiter, result->column_name_comp[i]);
  }
  fprintf(fp, "\n");

  for (i = 0; i < result->num_of_rows; i++) {
    fprintf(fp, "%.16g", *(value_time_p));
    value_time_p++;
    // Species
    for (j = 0; j < result->num_of_columns_sp; j++) {
      fprintf(fp, "%c%.16g", delimiter, *(value_sp_p));
      value_sp_p++;
    }
    // Parameters
    for (j = 0; j < result->num_of_columns_param; j++) {
      fprintf(fp, "%c%.16g", delimiter, *(value_param_p));
      value_param_p++;
    }
    // Compartments
    for (j = 0; j < result->num_of_columns_comp; j++) {
      fprintf(fp, "%c%.16g", delimiter, *(value_comp_p));
      value_comp_p++;
    }
    fprintf(fp, "\n");
  }
}

void write_separate_result(myResult* result, char* file_s, char* file_p, char* file_c) {
  FILE *fp_s, *fp_p, *fp_c;
  int i, j;
  char delimiter = ' ';
  double *value_time_p  = result->values_time;
  double *value_sp_p    = result->values_sp;
  double *value_param_p = result->values_param;
  double *value_comp_p  = result->values_comp;

  if ((fp_s = fopen(file_s, "w")) == NULL ) {
    fprintf(stderr, "Failed to open %s: %s\n", file_s, strerror(errno));
    return;
  }
  if ((fp_p = fopen(file_p, "w")) == NULL ) {
    fprintf(stderr, "Failed to open %s: %s\n", file_p, strerror(errno));
    return;
  }
  if ((fp_c = fopen(file_c, "w")) == NULL ) {
    fprintf(stderr, "Failed to open %s: %s\n", file_c, strerror(errno));
    return;
  }

  // Species
  for (i = 0; i < result->num_of_rows; i++) {
    fprintf(fp_s, "%.16g", *(value_time_p));
    value_time_p++;
    for (j = 0; j < result->num_of_columns_sp; j++) {
      fprintf(fp_s, "%c%.16g", delimiter, *(value_sp_p));
      value_sp_p++;
    }
    fprintf(fp_s, "\n");
  }
  // Parameters
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
  // Compartments
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
