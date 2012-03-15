/* libsbmlsim.i */
%module libsbmlsim
%{
#include "src/libsbmlsim/myResult.h"
  extern myResult* simulateSBMLFromString(const char *str, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method);
%}

/* %include "src/libsbmlsim/myResult.h" */
typedef _myResult myResult;
struct _myResult{
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
};

extern myResult* simulateSBMLFromString(const char *str, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method);
