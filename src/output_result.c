#include "sbml/common/common.h"
#include "sbml/SBMLReader.h"
#include "sbml/SBMLTypes.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include "header.h"

void print_result(myResult* result){
  output_result(result, stdout);
}

void print_result_to_file(myResult* result, char* file){
  FILE *fp;

  fp = fopen(file, "w");
  if (fp == NULL) {
    fprintf(stderr, "File %s open failed.\n", file);
  } else {
    output_result(result, fp);
    fclose(fp);
  }
}

void output_result(myResult* result, FILE* fp){
  int i, j;
  double *value_p = result->values;

  for (i = 0; i < result->num_of_columns; i++) {
    fprintf(fp, "%s ", result->column_name[i]);
  }
  fprintf(fp, "\n");
  for (i = 0; i < result->num_of_rows; i++) {
    for (j = 0; j < result->num_of_columns; j++) {
      fprintf(fp, "%.16g ", *(value_p));
      value_p++;
    }
    fprintf(fp, "\n");
  }
}
