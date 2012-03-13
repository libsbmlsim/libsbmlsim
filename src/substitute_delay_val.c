#include "libsbmlsim/libsbmlsim.h"

void substitute_delay_val(mySpecies *sp[], int num_of_species, myParameter *param[], int num_of_parameters, myCompartment *comp[], int num_of_compartments, myReaction *re[], int num_of_reactions, int cycle){
  int i, j;
  for(i=0; i<num_of_species; i++){
    if(sp[i]->delay_val != NULL){
      sp[i]->delay_val[cycle][0] = sp[i]->value;
      sp[i]->delay_val[cycle][1] = sp[i]->value;
      sp[i]->delay_val[cycle][2] = sp[i]->value;
      sp[i]->delay_val[cycle][3] = sp[i]->value;
    } 
  }
  for(i=0; i<num_of_parameters; i++){
    if(param[i]->delay_val != NULL){
      param[i]->delay_val[cycle][0] = param[i]->value;
      param[i]->delay_val[cycle][1] = param[i]->value;
      param[i]->delay_val[cycle][2] = param[i]->value;
      param[i]->delay_val[cycle][3] = param[i]->value;
    } 
  }
  for(i=0; i<num_of_compartments; i++){
    if(comp[i]->delay_val != NULL){
      comp[i]->delay_val[cycle][0] = comp[i]->value;
      comp[i]->delay_val[cycle][1] = comp[i]->value;
      comp[i]->delay_val[cycle][2] = comp[i]->value;
      comp[i]->delay_val[cycle][3] = comp[i]->value;
    } 
  }
  for(i=0; i<num_of_reactions; i++){
    for(j=0; j<re[i]->num_of_products; j++){
      if(re[i]->products[j]->delay_val != NULL){
	re[i]->products[j]->delay_val[cycle][0] = re[i]->products[j]->value;
	re[i]->products[j]->delay_val[cycle][1] = re[i]->products[j]->value;
	re[i]->products[j]->delay_val[cycle][2] = re[i]->products[j]->value;
	re[i]->products[j]->delay_val[cycle][3] = re[i]->products[j]->value;
      }
    }
    for(j=0; j<re[i]->num_of_reactants; j++){
      if(re[i]->reactants[j]->delay_val != NULL){
	re[i]->reactants[j]->delay_val[cycle][0] = re[i]->reactants[j]->value;
	re[i]->reactants[j]->delay_val[cycle][1] = re[i]->reactants[j]->value;
	re[i]->reactants[j]->delay_val[cycle][2] = re[i]->reactants[j]->value;
	re[i]->reactants[j]->delay_val[cycle][3] = re[i]->reactants[j]->value;
      }
    }
  }
}
