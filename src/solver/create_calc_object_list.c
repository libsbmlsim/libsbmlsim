/**
 * <!--------------------------------------------------------------------------
 * This file is part of libSBMLSim.  Please visit
 * http://fun.bio.keio.ac.jp/software/libsbmlsim/ for more
 * information about libSBMLSim and its latest version.
 *
 * Copyright (C) 2011-2013 by the Keio University, Yokohama, Japan
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution.
 * ---------------------------------------------------------------------- -->*/
#include "../libsbmlsim/libsbmlsim.h"

void create_calc_object_list(unsigned int num_of_species, unsigned int num_of_parameters, unsigned int num_of_compartments, unsigned int num_of_reactions, mySpecies *all_var_sp[], myParameter *all_var_param[], myCompartment *all_var_comp[], mySpeciesReference *all_var_spr[], mySpecies *var_sp[], myParameter *var_param[], myCompartment *var_comp[], mySpeciesReference *var_spr[], mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[]){

  unsigned int i, j, index1, index2;
  int is_included_in_reaction;

  index1 = 0;
  for(i=0; i<num_of_species; i++){
/*    if(!Species_getConstant(sp[i]->origin)){ */
      all_var_sp[index1++] = sp[i];
/*    } */
  }
  index1 = 0;
  for(i=0; i<num_of_parameters; i++){
    if(!Parameter_getConstant(param[i]->origin)){
      all_var_param[index1++] = param[i];
    }
  }
  index1 = 0;
  for(i=0; i<num_of_compartments; i++){
    if(!Compartment_getConstant(comp[i]->origin)){
      all_var_comp[index1++] = comp[i];
    }
  }
  index1 = 0;
  for(i=0; i<num_of_reactions; i++){
    for(j=0; j<re[i]->num_of_products; j++){
      if(!SpeciesReference_getConstant(re[i]->products[j]->origin)){
        all_var_spr[index1++] = re[i]->products[j];
      }
    }
    for(j=0; j<re[i]->num_of_reactants; j++){
      if(!SpeciesReference_getConstant(re[i]->reactants[j]->origin)){
        all_var_spr[index1++] = re[i]->reactants[j];
      }
    }
  }

  index2 = 0;
  for(i=0; i<num_of_species; i++){
    if(!Species_getConstant(sp[i]->origin)){
      if(sp[i]->depending_rule == NULL){
        is_included_in_reaction = 0;
        for(j=0; j<num_of_reactions; j++){
          if(!re[j]->is_fast){
            if(is_included_in_reaction){
              break;
            }
            for(index1=0; index1<re[j]->num_of_products; index1++){
              if(strcmp(Species_getId(sp[i]->origin), Species_getId(re[j]->products[index1]->mySp->origin)) == 0){
                is_included_in_reaction = 1;
              }
            }
            for(index1=0; index1<re[j]->num_of_reactants; index1++){
              if(strcmp(Species_getId(sp[i]->origin), Species_getId(re[j]->reactants[index1]->mySp->origin)) == 0){
                is_included_in_reaction = 1;
              }
            }
          }
        }
        if(is_included_in_reaction){
          var_sp[index2++] = sp[i];
        }
      }else if(sp[i]->depending_rule->is_rate){
        var_sp[index2++] = sp[i];
      }
    }
  }
  index1 = 0;
  for(i=0; i<num_of_parameters; i++){
    if(!Parameter_getConstant(param[i]->origin)){
      if(param[i]->depending_rule != NULL && param[i]->depending_rule->is_rate){
        var_param[index1++] = param[i];
      }
    }
  }
  index1 = 0;
  for(i=0; i<num_of_compartments; i++){
    if(!Compartment_getConstant(comp[i]->origin)){
      if(comp[i]->depending_rule != NULL && comp[i]->depending_rule->is_rate){
        var_comp[index1++] = comp[i];
      }
    }
  }
  index1 = 0;
  for(i=0; i<num_of_reactions; i++){
    for(j=0; j<re[i]->num_of_products; j++){
      if(!SpeciesReference_getConstant(re[i]->products[j]->origin)){
        if(re[i]->products[j]->depending_rule != NULL && re[i]->products[j]->depending_rule->is_rate){
          var_spr[index1++] = re[i]->products[j];
        }
      }
    }
    for(j=0; j<re[i]->num_of_reactants; j++){
      if(!SpeciesReference_getConstant(re[i]->reactants[j]->origin)){
        if(re[i]->reactants[j]->depending_rule != NULL && re[i]->reactants[j]->depending_rule->is_rate){
          var_spr[index1++] = re[i]->reactants[j];
        }
      }
    }
  }

}
