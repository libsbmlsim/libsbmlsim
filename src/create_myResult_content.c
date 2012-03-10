#include "sbml/common/common.h"
#include "sbml/SBMLReader.h"
#include "sbml/SBMLTypes.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include "header.h"

//create contents of myResult object
void create_myResult_content(Model_t *m, myResult* result, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], double sim_time, double dt, int print_interval) {
  int i, last;
  int num_of_species = Model_getNumSpecies(m);
  int num_of_parameters = Model_getNumParameters(m);
  int num_of_compartments = Model_getNumCompartments(m);
  int end_cycle = get_end_cycle(sim_time, dt);

  //XXX
  result->num_of_columns = num_of_species + num_of_parameters + num_of_compartments + 1;
  result->num_of_columns_sp = num_of_species;
  result->num_of_columns_param = num_of_parameters;
  result->num_of_columns_comp = num_of_compartments;
  result->num_of_rows = end_cycle / print_interval + 1;
  //XXX
  result->column_name = (const char **)malloc(sizeof(char *) * result->num_of_columns);
  result->column_name_time  = "#time";
  result->column_name_sp    = (const char **)malloc(sizeof(char *) * num_of_species);
  result->column_name_param = (const char **)malloc(sizeof(char *) * num_of_parameters);
  result->column_name_comp  = (const char **)malloc(sizeof(char *) * num_of_compartments);
  //XXX
  result->values = (double *)malloc(sizeof(double) * result->num_of_columns * result->num_of_rows);
  result->values_time = (double *)malloc(sizeof(double) * result->num_of_rows);
  result->values_sp = (double *)malloc(sizeof(double) * num_of_species * result->num_of_rows);
  result->values_param = (double *)malloc(sizeof(double) * num_of_parameters * result->num_of_rows);
  result->values_comp = (double *)malloc(sizeof(double) * num_of_compartments * result->num_of_rows);
  result->column_name[0] = "#time";
  last = 1;
  for(i=0; i<num_of_species; i++){
    result->column_name[last] = Species_getId(mySp[i]->origin); // XXX
    result->column_name_sp[i] = Species_getId(mySp[i]->origin);
    last++;
  }
  for(i=0; i<num_of_parameters; i++){
    result->column_name[last] = Parameter_getId(myParam[i]->origin);
    result->column_name_param[i] = Parameter_getId(myParam[i]->origin);
    last++;
  }
  for(i=0; i<num_of_compartments; i++){
    result->column_name[last] = Compartment_getId(myComp[i]->origin);
    result->column_name_comp[i] = Compartment_getId(myComp[i]->origin);
    last++;
  }
  return;
}
