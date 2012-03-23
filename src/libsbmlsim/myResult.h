#ifndef LibSBMLSim_MyResult_h
#define LibSBMLSim_MyResult_h

#include "errorcodes.h"

typedef struct _myResult myResult;

struct _myResult{
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
};

#endif /* LibSBMLSim_MyResult_h */
