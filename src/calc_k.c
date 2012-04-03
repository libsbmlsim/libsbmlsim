/**
 * <!--------------------------------------------------------------------------
 * This file is part of libSBMLSim.  Please visit
 * http://fun.bio.keio.ac.jp/software/libsbmlsim/ for more
 * information about libSBMLSim and its latest version.
 *
 * Copyright (C) 2011-2012 by the Keio University, Yokohama, Japan
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution.
 * ---------------------------------------------------------------------- -->*/
#include "libsbmlsim/libsbmlsim.h"

void calc_k(mySpecies *sp[], unsigned int sp_num, myParameter *param[], unsigned int param_num, myCompartment *comp[], unsigned int comp_num, mySpeciesReference *spr[], unsigned int spr_num, myReaction *re[], unsigned int re_num, myRule *rule[], unsigned int rule_num, int cycle, double dt, double *reverse_time, int use_rk, int call_first_time_in_cycle){
  unsigned int i, j;
  double k = 0;
  double rk_cef[4] = {0.5, 0.5, 1, 0};
  int step;
  int step_num;

  if(use_rk){
    step_num = 4;
  }else{
    step_num = 1;
  }

  for(step=0; step<step_num; step++){
    if(call_first_time_in_cycle){
      /* substitute to delay buf */
      for(i=0; i<sp_num; i++){
        if(sp[i]->delay_val != NULL){
          sp[i]->delay_val[cycle][step] = sp[i]->temp_value;
        }
      }
      for(i=0; i<param_num; i++){
        if(param[i]->delay_val != NULL){
          param[i]->delay_val[cycle][step] = param[i]->temp_value;
        }
      }
      for(i=0; i<comp_num; i++){
        if(comp[i]->delay_val != NULL){
          comp[i]->delay_val[cycle][step] = comp[i]->temp_value;
        }
      }
      for(i=0; i<spr_num; i++){
        if(spr[i]->delay_val != NULL){
          spr[i]->delay_val[cycle][step] = spr[i]->temp_value;
        }
      }
    }
    /* initialize k */
    for(i=0; i<sp_num; i++){
      sp[i]->k[step] = 0;
    }
    for(i=0; i<param_num; i++){
      param[i]->k[step] = 0;
    }
    for(i=0; i<comp_num; i++){
      comp[i]->k[step] = 0;
    }
    for(i=0; i<spr_num; i++){
      spr[i]->k[step] = 0;
    }
    /* reaction */
    for(i=0; i<re_num; i++){
      if(!re[i]->is_fast){
        k = calc(re[i]->eq, dt, cycle, reverse_time, step);
        for(j=0; j<re[i]->num_of_products; j++){
          if(!Species_getBoundaryCondition(re[i]->products[j]->mySp->origin)){
            re[i]->products[j]->mySp->k[step] += calc(re[i]->products[j]->eq, dt, cycle, reverse_time, step)*k; 
          }
        }
        for(j=0; j<re[i]->num_of_reactants; j++){
          if(!Species_getBoundaryCondition(re[i]->reactants[j]->mySp->origin)){
            re[i]->reactants[j]->mySp->k[step] -= calc(re[i]->reactants[j]->eq, dt, cycle, reverse_time, step)*k; 
          }
        }
      }
    }
    /* rule */
    for(i=0; i<rule_num; i++){
      if(rule[i]->target_species != NULL){/*  calculate math in rule */
        rule[i]->target_species->k[step] += calc(rule[i]->eq, dt, cycle, reverse_time, step);
      }else if(rule[i]->target_parameter != NULL){/*  calculate math in rule */
        rule[i]->target_parameter->k[step] += calc(rule[i]->eq, dt, cycle, reverse_time, step);
      }else if(rule[i]->target_compartment != NULL){/*  calculate math in rule */
        rule[i]->target_compartment->k[step] += calc(rule[i]->eq, dt, cycle, reverse_time, step);
      }else if(rule[i]->target_species_reference != NULL){/*  calculate math in rule */
        rule[i]->target_species_reference->k[step] += calc(rule[i]->eq, dt, cycle, reverse_time, step);
      }
    }

    if(use_rk){
      /* species */
      for(i=0; i<sp_num; i++){
        if(sp[i]->depending_rule != NULL && sp[i]->depending_rule->is_assignment){
          sp[i]->temp_value = sp[i]->k[step];
        }else{
          sp[i]->temp_value = sp[i]->value + sp[i]->k[step]*dt*rk_cef[step];
        }
      }
      /* parameter */
      for(i=0; i<param_num; i++){
        if(param[i]->depending_rule != NULL && param[i]->depending_rule->is_assignment){
          param[i]->temp_value = param[i]->k[step];
        }else{
          param[i]->temp_value = param[i]->value + param[i]->k[step]*dt*rk_cef[step];
        }
      }
      /* compartment */
      for(i=0; i<comp_num; i++){
        if(comp[i]->depending_rule != NULL && comp[i]->depending_rule->is_assignment){
          /* new code */
          for(j=0; j<comp[i]->num_of_including_species; j++){
            if(comp[i]->including_species[j]->is_concentration){
              comp[i]->including_species[j]->temp_value = comp[i]->including_species[j]->temp_value*comp[i]->temp_value/comp[i]->k[step];
            }
          }
         /* new code end */
          comp[i]->temp_value = comp[i]->k[step];
        }else{
          /* new code */
          for(j=0; j<comp[i]->num_of_including_species; j++){
            if(comp[i]->including_species[j]->is_concentration){
              comp[i]->including_species[j]->temp_value = comp[i]->including_species[j]->temp_value*comp[i]->temp_value/(comp[i]->value + comp[i]->k[step]*dt*rk_cef[step]);
            }
          }
         /* new code end */
          comp[i]->temp_value = comp[i]->value + comp[i]->k[step]*dt*rk_cef[step];
        }
      }
      /* species reference */
      for(i=0; i<spr_num; i++){
        if(spr[i]->depending_rule != NULL && spr[i]->depending_rule->is_assignment){
          spr[i]->temp_value = spr[i]->k[step];
        }else{
          spr[i]->temp_value = spr[i]->value + spr[i]->k[step]*dt*rk_cef[step];
        }
      }      
    } 
  }
}
