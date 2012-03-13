#include "libsbmlsim/libsbmlsim.h"

void calc_temp_value(mySpecies *sp[], int sp_num, myParameter *param[], int param_num, myCompartment *comp[], int comp_num, mySpeciesReference *spr[], int spr_num, int cycle, double dt, int use_rk){
  int i, j;
  
  if(use_rk){
    //species
    for(i=0; i<sp_num; i++){
      if(sp[i]->depending_rule != NULL && !sp[i]->depending_rule->is_rate){
	sp[i]->temp_value = (sp[i]->k[0]+2*sp[i]->k[1]+2*sp[i]->k[2]+sp[i]->k[3])/6;
      }else{
	sp[i]->temp_value = sp[i]->value + (sp[i]->k[0]+2*sp[i]->k[1]+2*sp[i]->k[2]+sp[i]->k[3])/6*dt;
      }
    }
    //parameter
    for(i=0; i<param_num; i++){
      if(param[i]->depending_rule != NULL && !param[i]->depending_rule->is_rate){
	param[i]->temp_value = (param[i]->k[0]+2*param[i]->k[1]+2*param[i]->k[2]+param[i]->k[3])/6;
      }else{
	param[i]->temp_value = param[i]->value + (param[i]->k[0]+2*param[i]->k[1]+2*param[i]->k[2]+param[i]->k[3])/6*dt;
      }
    }
    //compartment
    for(i=0; i<comp_num; i++){
      if(comp[i]->depending_rule != NULL && !comp[i]->depending_rule->is_rate){
	comp[i]->temp_value = (comp[i]->k[0]+2*comp[i]->k[1]+2*comp[i]->k[2]+comp[i]->k[3])/6;
      }else{
	comp[i]->temp_value = comp[i]->value + (comp[i]->k[0]+2*comp[i]->k[1]+2*comp[i]->k[2]+comp[i]->k[3])/6*dt;
      }
      //new code
      for(j=0; j<comp[i]->num_of_including_species; j++){
        if(comp[i]->including_species[j]->is_concentration){
          comp[i]->including_species[j]->temp_value = comp[i]->including_species[j]->temp_value*comp[i]->value/comp[i]->temp_value;
        }
      }
      //
    }
    //species reference
    for(i=0; i<spr_num; i++){
      if(spr[i]->depending_rule != NULL && !spr[i]->depending_rule->is_rate){
	spr[i]->temp_value = (spr[i]->k[0]+2*spr[i]->k[1]+2*spr[i]->k[2]+spr[i]->k[3])/6;
      }else{
	spr[i]->temp_value = spr[i]->value + (spr[i]->k[0]+2*spr[i]->k[1]+2*spr[i]->k[2]+spr[i]->k[3])/6*dt;
      }
    }
  }else{
    //species
    for(i=0; i<sp_num; i++){
      if(sp[i]->depending_rule != NULL && !sp[i]->depending_rule->is_rate){
	sp[i]->temp_value = sp[i]->k[0];
      }else{
	sp[i]->temp_value = sp[i]->value + sp[i]->k[0]*dt;
      }
    }
    //parameter
    for(i=0; i<param_num; i++){
      if(param[i]->depending_rule != NULL && !param[i]->depending_rule->is_rate){
	param[i]->temp_value = param[i]->k[0];
      }else{
	param[i]->temp_value = param[i]->value + param[i]->k[0]*dt;
      }
    }
    //compartment
    for(i=0; i<comp_num; i++){
      if(comp[i]->depending_rule != NULL && !comp[i]->depending_rule->is_rate){
	comp[i]->temp_value = comp[i]->k[0];
      }else{
	comp[i]->temp_value = comp[i]->value + comp[i]->k[0]*dt;
      }
      //new code
      for(j=0; j<comp[i]->num_of_including_species; j++){
        if(comp[i]->including_species[j]->is_concentration){
          comp[i]->including_species[j]->temp_value = comp[i]->including_species[j]->temp_value*comp[i]->value/comp[i]->temp_value;
        }
      }
      //
    }
    //species reference
    for(i=0; i<spr_num; i++){
      if(spr[i]->depending_rule != NULL && !spr[i]->depending_rule->is_rate){
	spr[i]->temp_value = spr[i]->k[0];
      }else{
	spr[i]->temp_value = spr[i]->value + spr[i]->k[0]*dt;
      }
    }
  }
  
}
