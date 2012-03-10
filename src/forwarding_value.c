#include "sbml/common/common.h"
#include "sbml/SBMLReader.h"
#include "sbml/SBMLTypes.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include "header.h"

void forwarding_value(mySpecies *sp[], int sp_num, myParameter *param[], int param_num, myCompartment *comp[], int comp_num, mySpeciesReference *spr[], int spr_num){
  int i;
  for(i=0; i<sp_num; i++){
    sp[i]->value = sp[i]->temp_value;
  }
  for(i=0; i<param_num; i++){
    param[i]->value = param[i]->temp_value;
  }
  for(i=0; i<comp_num; i++){
    comp[i]->value = comp[i]->temp_value;
  }
  for(i=0; i<spr_num; i++){
    spr[i]->value = spr[i]->temp_value;
  }  
}
