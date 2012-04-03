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

void calc_temp_value(mySpecies *sp[], unsigned int sp_num, myParameter *param[], unsigned int param_num, myCompartment *comp[], unsigned int comp_num, mySpeciesReference *spr[], unsigned int spr_num, double dt, int use_rk){
  unsigned int i, j;

  if(use_rk){
    /* species */
    for(i=0; i<sp_num; i++){
      if(sp[i]->depending_rule != NULL && !sp[i]->depending_rule->is_rate){
        sp[i]->temp_value = (sp[i]->k[0]+2*sp[i]->k[1]+2*sp[i]->k[2]+sp[i]->k[3])/6;
      }else{
        sp[i]->temp_value = sp[i]->value + (sp[i]->k[0]+2*sp[i]->k[1]+2*sp[i]->k[2]+sp[i]->k[3])/6*dt;
      }
    }
    /* parameter */
    for(i=0; i<param_num; i++){
      if(param[i]->depending_rule != NULL && !param[i]->depending_rule->is_rate){
        param[i]->temp_value = (param[i]->k[0]+2*param[i]->k[1]+2*param[i]->k[2]+param[i]->k[3])/6;
      }else{
        param[i]->temp_value = param[i]->value + (param[i]->k[0]+2*param[i]->k[1]+2*param[i]->k[2]+param[i]->k[3])/6*dt;
      }
    }
    /* compartment */
    for(i=0; i<comp_num; i++){
      if(comp[i]->depending_rule != NULL && !comp[i]->depending_rule->is_rate){
        comp[i]->temp_value = (comp[i]->k[0]+2*comp[i]->k[1]+2*comp[i]->k[2]+comp[i]->k[3])/6;
      }else{
        comp[i]->temp_value = comp[i]->value + (comp[i]->k[0]+2*comp[i]->k[1]+2*comp[i]->k[2]+comp[i]->k[3])/6*dt;
      }
      /* new code */
      for(j=0; j<comp[i]->num_of_including_species; j++){
        if(comp[i]->including_species[j]->is_concentration){
          comp[i]->including_species[j]->temp_value = comp[i]->including_species[j]->temp_value*comp[i]->value/comp[i]->temp_value;
        }
      }
     /* new code end */
    }
    /* species reference */
    for(i=0; i<spr_num; i++){
      if(spr[i]->depending_rule != NULL && !spr[i]->depending_rule->is_rate){
        spr[i]->temp_value = (spr[i]->k[0]+2*spr[i]->k[1]+2*spr[i]->k[2]+spr[i]->k[3])/6;
      }else{
        spr[i]->temp_value = spr[i]->value + (spr[i]->k[0]+2*spr[i]->k[1]+2*spr[i]->k[2]+spr[i]->k[3])/6*dt;
      }
    }
  }else{
    /* species */
    for(i=0; i<sp_num; i++){
      if(sp[i]->depending_rule != NULL && !sp[i]->depending_rule->is_rate){
        sp[i]->temp_value = sp[i]->k[0];
      }else{
        sp[i]->temp_value = sp[i]->value + sp[i]->k[0]*dt;
      }
    }
    /* parameter */
    for(i=0; i<param_num; i++){
      if(param[i]->depending_rule != NULL && !param[i]->depending_rule->is_rate){
        param[i]->temp_value = param[i]->k[0];
      }else{
        param[i]->temp_value = param[i]->value + param[i]->k[0]*dt;
      }
    }
    /* compartment */
    for(i=0; i<comp_num; i++){
      if(comp[i]->depending_rule != NULL && !comp[i]->depending_rule->is_rate){
        comp[i]->temp_value = comp[i]->k[0];
      }else{
        comp[i]->temp_value = comp[i]->value + comp[i]->k[0]*dt;
      }
      /* new code */
      for(j=0; j<comp[i]->num_of_including_species; j++){
        if(comp[i]->including_species[j]->is_concentration){
          comp[i]->including_species[j]->temp_value = comp[i]->including_species[j]->temp_value*comp[i]->value/comp[i]->temp_value;
        }
      }
     /* new code end */
    }
    /* species reference */
    for(i=0; i<spr_num; i++){
      if(spr[i]->depending_rule != NULL && !spr[i]->depending_rule->is_rate){
        spr[i]->temp_value = spr[i]->k[0];
      }else{
        spr[i]->temp_value = spr[i]->value + spr[i]->k[0]*dt;
      }
    }
  }

}
