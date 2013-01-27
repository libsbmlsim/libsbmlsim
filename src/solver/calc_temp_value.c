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

double calc_eps(double value) {
  if (1 + value == value) {
    value *= 2.0;
    calc_eps(value);
	}
	return value;
}

double calc_error(double dxdt, double dxdt4, double cur_value, double next_value, double atol, double rtol) {
  double sci;
  double err;

  sci = atol + my_fmax(fabs(cur_value), fabs(next_value)) * rtol;
  err = (dxdt4 - dxdt) / sci;
  return err * err;
}

double calc_sum_error(mySpecies *sp[], unsigned int sp_num, myParameter *param[], unsigned int param_num, myCompartment *comp[], unsigned int comp_num, mySpeciesReference *spr[], unsigned int spr_num, double dt, int use_rk, double atol, double rtol, int* ode_num, int time_progressed_flag, int order){
  unsigned int i, j;
  int l,m;
  int step;
  int step_num = 6;
  double sum_error = 0.0;
  double rk_ce[2][6];

  /* Runge-Kutta-Fehlberg table */
  double rk_ce_f[2][6] = {{16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0},
						  {25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0}};
  /* Cash-Karp table */
  double rk_ce_c[2][6] = {{2825.0/27648.0, 0.0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0, 1.0/4.0},
						  {37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0}};

  /* allocate memory */
  double *sp_dxdt  = (double *)malloc(sizeof(double) * sp_num);
  double *sp_dxdt4 = (double *)malloc(sizeof(double) * sp_num);
  double *param_dxdt  = (double *)malloc(sizeof(double) * param_num);
  double *param_dxdt4 = (double *)malloc(sizeof(double) * param_num);
  double *comp_dxdt  = (double *)malloc(sizeof(double) * comp_num);
  double *comp_dxdt4 = (double *)malloc(sizeof(double) * comp_num);
  double *spr_dxdt  = (double *)malloc(sizeof(double) * spr_num);
  double *spr_dxdt4 = (double *)malloc(sizeof(double) * spr_num);

  /* prepare Butcher tableau */
  if(order == 5){
	  for(l=0; l<2; l++){
		  for(m=0; m<step_num; m++){
			  rk_ce[l][m] = rk_ce_f[l][m];
		  }
	  }
  }else if (order == 6){
	  for(l=0; l<2; l++){
		  for(m=0; m<step_num; m++){
			  rk_ce[l][m] = rk_ce_c[l][m];
		  }
	  }
  }else{
	  fprintf(stderr, "Butcher tableau is implicit.\n");
	  exit(1);
  }

  if(use_rk){
	  /* species */
	  for(i=0; i<sp_num; i++){
		  sp_dxdt[i] = 0.0;
		  sp_dxdt4[i] = 0.0;
		  for(step=0; step<step_num; step++){
			  sp_dxdt[i] += rk_ce[0][step]*sp[i]->k[step];
			  sp_dxdt4[i] += rk_ce[1][step]*sp[i]->k[step];
		  }
		  sum_error += calc_error(sp_dxdt[i]*dt, sp_dxdt4[i]*dt, sp[i]->value, sp[i]->value + sp_dxdt4[i]*dt, atol, rtol);
	  }
	  /* parameter */
	  for(i=0; i<param_num; i++){
		  param_dxdt[i] = 0.0;
		  param_dxdt4[i] = 0.0;
		  for(step=0; step<step_num; step++){
			  param_dxdt[i] += rk_ce[0][step]*param[i]->k[step];
			  param_dxdt4[i] += rk_ce[1][step]*param[i]->k[step];
		  }
		  sum_error += calc_error(param_dxdt[i]*dt, param_dxdt4[i]*dt, param[i]->value, param[i]->value + param_dxdt4[i]*dt, atol, rtol);
	  }
	  /* compartment */
	  for(i=0; i<comp_num; i++){
		  comp_dxdt[i] = 0.0;
		  comp_dxdt4[i] = 0.0;
		  for(step=0; step<step_num; step++){
			  comp_dxdt[i] += rk_ce[0][step]*comp[i]->k[step];
			  comp_dxdt4[i] += rk_ce[1][step]*comp[i]->k[step];
		  }
		  sum_error += calc_error(comp_dxdt[i]*dt, comp_dxdt4[i]*dt, comp[i]->value, comp[i]->value + comp_dxdt4[i]*dt, atol, rtol);
	  }
	  /* species reference */
	  for(i=0; i<spr_num; i++){
		  spr_dxdt[i] = 0.0;
		  spr_dxdt4[i] = 0.0;
		  for(step=0; step<step_num; step++){
			  spr_dxdt[i] += rk_ce[0][step]*spr[i]->k[step];
			  spr_dxdt4[i] += rk_ce[1][step]*spr[i]->k[step];
		  }
		  sum_error += calc_error(spr_dxdt[i]*dt, spr_dxdt4[i]*dt, spr[i]->value, spr[i]->value + spr_dxdt4[i]*dt, atol, rtol);
	  }

	  if (*(ode_num) == 0){
      /* Don't forget to free allocated memory before return! */
      free_dxdt(sp_dxdt, sp_dxdt4, param_dxdt, param_dxdt4, comp_dxdt, comp_dxdt4, spr_dxdt, spr_dxdt4);
      return 0;
	  }else {
		  sum_error = sqrt(sum_error / *(ode_num));
	  }

	  if (sum_error > 1.0){
		  for(i=0; i<sp_num; i++){ /* reset */
			  sp[i]->temp_value = sp[i]->value;
		  }
		  for(i=0; i<param_num; i++){ /* reset */
			  param[i]->temp_value = param[i]->value;
		  }
		  for(i=0; i<comp_num; i++){ /* reset */
			  comp[i]->temp_value = comp[i]->value;
		  }
		  for(i=0; i<spr_num; i++){ /* reset */
			  spr[i]->temp_value = spr[i]->value;
		  }
	  }

	  if(time_progressed_flag){
		  /* species */
		  for(i=0; i<sp_num; i++){
			  if(sp[i]->depending_rule != NULL && !sp[i]->depending_rule->is_rate){
				  sp[i]->temp_value = sp_dxdt[i];
			  }else{
				  sp[i]->temp_value = sp[i]->value + sp_dxdt[i]*dt;
			  }
		  }
		  /* parameter */
		  for(i=0; i<param_num; i++){
			  if(param[i]->depending_rule != NULL && !param[i]->depending_rule->is_rate){
				  param[i]->temp_value = param_dxdt[i];
			  }else{
				  param[i]->temp_value = param[i]->value + param_dxdt[i]*dt;
			  }
		  }
		  /* compartment */
		  for(i=0; i<comp_num; i++){
			  if(comp[i]->depending_rule != NULL && !comp[i]->depending_rule->is_rate){
				  comp[i]->temp_value = comp_dxdt[i];
			  }else{
				  comp[i]->temp_value = comp[i]->value + comp_dxdt[i]*dt;
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
				  spr[i]->temp_value = spr_dxdt[i];
			  }else{
				  spr[i]->temp_value = spr[i]->value + spr_dxdt[i]*dt;
			  }
		  }
	  }
  }
  /* Don't forget to free allocated memory before return! */
  free_dxdt(sp_dxdt, sp_dxdt4, param_dxdt, param_dxdt4, comp_dxdt, comp_dxdt4, spr_dxdt, spr_dxdt4);
  return sum_error;
}

void free_dxdt(double *sp_dxdt, double *sp_dxdt4, double *param_dxdt, double *param_dxdt4, double *comp_dxdt, double *comp_dxdt4, double *spr_dxdt, double *spr_dxdt4) { 
  free(sp_dxdt);
  free(sp_dxdt4);
  free(param_dxdt);
  free(param_dxdt4);
  free(comp_dxdt);
  free(comp_dxdt4);
  free(spr_dxdt);
  free(spr_dxdt4);
}
