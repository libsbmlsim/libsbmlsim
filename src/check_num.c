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
#include "libsbmlsim/libsbmlsim.h"

void check_num(unsigned int num_of_species, unsigned int num_of_parameters, unsigned int num_of_compartments, unsigned int num_of_reactions, unsigned int *num_of_all_var_species, unsigned int *num_of_all_var_parameters, unsigned int *num_of_all_var_compartments, unsigned int *num_of_all_var_species_reference, unsigned int *num_of_var_species, unsigned int *num_of_var_parameters, unsigned int *num_of_var_compartments, unsigned int *num_of_var_species_reference, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[]){

  unsigned int i, j, index;
  int is_included_in_reaction;

  for(i=0; i<num_of_species; i++){
/*    if(!Species_getConstant(sp[i]->origin)){ */
      (*num_of_all_var_species)++;
/*    } */
  }
  for(i=0; i<num_of_parameters; i++){
    if(!Parameter_getConstant(param[i]->origin)){
      (*num_of_all_var_parameters)++;
    }
  }
  for(i=0; i<num_of_compartments; i++){
    if(!Compartment_getConstant(comp[i]->origin)){
      (*num_of_all_var_compartments)++;
    }
  }
  for(i=0; i<num_of_reactions; i++){
    for(j=0; j<re[i]->num_of_products; j++){
      if(!SpeciesReference_getConstant(re[i]->products[j]->origin)){
        (*num_of_all_var_species_reference)++;
      }
    }
    for(j=0; j<re[i]->num_of_reactants; j++){
      if(!SpeciesReference_getConstant(re[i]->reactants[j]->origin)){
        (*num_of_all_var_species_reference)++;
      }
    }
  }

  for(i=0; i<num_of_species; i++){
    if(!Species_getConstant(sp[i]->origin)){
      if(sp[i]->depending_rule == NULL){
        is_included_in_reaction = 0;
        for(j=0; j<num_of_reactions; j++){
          if(!re[j]->is_fast){
            if(is_included_in_reaction){
              break;
            }
            for(index=0; index<re[j]->num_of_products; index++){
              if(strcmp(Species_getId(sp[i]->origin), Species_getId(re[j]->products[index]->mySp->origin)) == 0){
                is_included_in_reaction = 1;
              }
            }
            for(index=0; index<re[j]->num_of_reactants; index++){
              if(strcmp(Species_getId(sp[i]->origin), Species_getId(re[j]->reactants[index]->mySp->origin)) == 0){
                is_included_in_reaction = 1;
              }
            }
          }
        }
        if(is_included_in_reaction){
          (*num_of_var_species)++;
        }
      }else if(sp[i]->depending_rule->is_rate){
        (*num_of_var_species)++;
      }
    }
  }
  for(i=0; i<num_of_parameters; i++){
    if(!Parameter_getConstant(param[i]->origin)){
      if(param[i]->depending_rule != NULL && param[i]->depending_rule->is_rate){
        (*num_of_var_parameters)++;
      }
    }
  }
  for(i=0; i<num_of_compartments; i++){
    if(!Compartment_getConstant(comp[i]->origin)){
      if(comp[i]->depending_rule != NULL && comp[i]->depending_rule->is_rate){
        (*num_of_var_compartments)++;
      }
    }
  }
  for(i=0; i<num_of_reactions; i++){
    for(j=0; j<re[i]->num_of_products; j++){
      if(!SpeciesReference_getConstant(re[i]->products[j]->origin)){
        if(re[i]->products[j]->depending_rule != NULL && re[i]->products[j]->depending_rule->is_rate){
          (*num_of_var_species_reference)++;
        }
      }
    }
    for(j=0; j<re[i]->num_of_reactants; j++){
      if(!SpeciesReference_getConstant(re[i]->reactants[j]->origin)){
        if(re[i]->reactants[j]->depending_rule != NULL && re[i]->reactants[j]->depending_rule->is_rate){
          (*num_of_var_species_reference)++;
        }
      }
    }
  }

}
