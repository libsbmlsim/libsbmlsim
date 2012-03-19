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
    /*result->column_name_sp[i] = dupstr(Species_getId(mySp[i]->origin));*/
    result->column_name_sp[i] = Species_getId(mySp[i]->origin);
  }
  for(i=0; i<num_of_parameters; i++){
    /*result->column_name_param[i] = dupstr(Parameter_getId(myParam[i]->origin));*/
    result->column_name_param[i] = Parameter_getId(myParam[i]->origin);
  }
  for(i=0; i<num_of_compartments; i++){
    /*result->column_name_comp[i] = dupstr(Compartment_getId(myComp[i]->origin));*/
    result->column_name_comp[i] = Compartment_getId(myComp[i]->origin);
  }
  return result;
}

void free_myResult(myResult *res)
{
  return;
  int i;
  free((void *)res->column_name_time);
  for (i = 0; i < res->num_of_columns_sp; i++)
    free((void *)res->column_name_sp[i]);
  free(res->column_name_sp);
  for (i = 0; i < res->num_of_columns_param; i++)
    free((void *)res->column_name_param[i]);
  free(res->column_name_param);
  for (i = 0; i < res->num_of_columns_comp; i++)
    free((void *)res->column_name_comp[i]);
  free(res->column_name_comp);
  free(res->values_time);
  free(res->values_sp);
  free(res->values_param);
  free(res->values_comp);
  free(res);
}
