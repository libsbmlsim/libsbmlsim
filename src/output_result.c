#include "sbml/common/common.h"
#include "sbml/SBMLReader.h"
#include "sbml/SBMLTypes.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include "header.h"

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

void output_result(myResult* result, FILE* fp, char delimiter){
  int i, j;
  double *value_p = result->values;

  for (i = 0; i < result->num_of_columns; i++) {
    fprintf(fp, "%s%c", result->column_name[i], delimiter);
  }
  fprintf(fp, "\n");
  for (i = 0; i < result->num_of_rows; i++) {
    for (j = 0; j < result->num_of_columns; j++) {
      fprintf(fp, "%.16g%c", *(value_p), delimiter);
      value_p++;
    }
    fprintf(fp, "\n");
  }
}
