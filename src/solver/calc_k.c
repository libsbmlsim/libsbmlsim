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

void calc_kf(mySpecies *sp[], unsigned int sp_num, myParameter *param[], unsigned int param_num, myCompartment *comp[], unsigned int comp_num, mySpeciesReference *spr[], unsigned int spr_num, myReaction *re[], unsigned int re_num, myRule *rule[], unsigned int rule_num, int cycle, double dt, double *reverse_time, int use_rk, int call_first_time_in_cycle, double* time, myResult* res, myAlgebraicEquations *algEq, int print_interval, int* err_zero_flag, int order){
  unsigned int i, j;
  int l, m;
  double k = 0;
  double rk_ce[5][5];
  double rk_cet[6];
  int step;
  int step_num;
  int error;
  double time_step[6];

  /* Runge-Kutta-Fehlberg table */
  double rk_ce_f[5][5] = {{1.0/4.0, 0.0, 0.0, 0.0, 0.0},
						  {3.0/32.0, 9.0/32.0, 0.0, 0.0, 0.0},
						  {1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0, 0.0},
						  {439.0/216.0, -8.0,3680.0/513.0, -845.0/4104.0, 0.0},
						  {-8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0}};
  double rk_cet_f[6] = {0.0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0};
  /* Cash-Karp table */
  double rk_ce_c[5][5] = {{1.0/4.0, 0.0, 0.0, 0.0, 0.0},
						  {3.0/32.0, 9.0/32.0, 0.0, 0.0, 0.0},
						  {1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0, 0.0},
						  {439.0/216.0, -8.0,3680.0/513.0, -845.0/4104.0, 0.0},
						  {-8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0}};
  double rk_cet_c[6] = {0.0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0};

  if(use_rk){
    step_num = 6;
  }else{
    step_num = 1;
  }

  /* prepare Butcher tableau */
  if (order == 5){
	  for(l=0; l<step_num-1; l++){
		  for(m=0; m<step_num-1; m++){
			  rk_ce[l][m] = rk_ce_f[l][m];
		  }
	  }
	  for(l=0; l<step_num; l++){
		  rk_cet[l] = rk_cet_f[l];
	  }
  }else if(order == 6){
	  for(l=0; l<step_num-1; l++){
		  for(m=0; m<step_num-1; m++){
			  rk_ce[l][m] = rk_ce_c[l][m];
		  }
	  }
	  for(l=0; l<step_num; l++){
		  rk_cet[l] = rk_cet_c[l];
	  }
  }else{
	  fprintf(stderr, "Butcher tableau is implicit.\n");
	  exit(1);
  }

  for(step=0; step<step_num; step++){
	  time_step[step] = *time + rk_cet[step]*dt;
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
		  k = calcf(re[i]->eq, dt, cycle, reverse_time, step, time, &time_step[step], res, print_interval, err_zero_flag);
        for(j=0; j<re[i]->num_of_products; j++){
          if(!Species_getBoundaryCondition(re[i]->products[j]->mySp->origin)){
			  re[i]->products[j]->mySp->k[step] += calcf(re[i]->products[j]->eq, dt, cycle, reverse_time, step, time, &time_step[step], res, print_interval, err_zero_flag)*k;
          }
        }
        for(j=0; j<re[i]->num_of_reactants; j++){
          if(!Species_getBoundaryCondition(re[i]->reactants[j]->mySp->origin)){
			  re[i]->reactants[j]->mySp->k[step] -= calcf(re[i]->reactants[j]->eq, dt, cycle, reverse_time, step, time, &time_step[step], res, print_interval, err_zero_flag)*k;
          }
        }
      }
    }

    /* rule */
    for(i=0; i<rule_num; i++){
      if(rule[i]->target_species != NULL){/*  calculate math in rule */
 		  rule[i]->target_species->k[step] += calcf(rule[i]->eq, dt, cycle, reverse_time, step, time, &time_step[step], res, print_interval, err_zero_flag);
      }else if(rule[i]->target_parameter != NULL){/*  calculate math in rule */
 		  rule[i]->target_parameter->k[step] += calcf(rule[i]->eq, dt, cycle, reverse_time, step, time, &time_step[step], res, print_interval, err_zero_flag);
      }else if(rule[i]->target_compartment != NULL){/*  calculate math in rule */
 		  rule[i]->target_compartment->k[step] += calcf(rule[i]->eq, dt, cycle, reverse_time, step, time, &time_step[step], res, print_interval, err_zero_flag);
      }else if(rule[i]->target_species_reference != NULL){/*  calculate math in rule */
 		  rule[i]->target_species_reference->k[step] += calcf(rule[i]->eq, dt, cycle, reverse_time, step, time, &time_step[step], res, print_interval, err_zero_flag);
      }
    }


	if(use_rk){
		/* species */
		for(i=0; i<sp_num; i++){
			if(sp[i]->depending_rule != NULL && sp[i]->depending_rule->is_assignment){
				sp[i]->temp_value = sp[i]->k[step];
			}else{
				sp[i]->temp_value = sp[i]->value;
				for(l=0; l<step; l++){
					sp[i]->temp_value += sp[i]->k[l]*dt*rk_ce[step-1][l];
				}
				error = calc_by_algebraic(algEq, cycle, dt, *reverse_time, time, res, print_interval, err_zero_flag);
				if (error == 1) {
					fprintf(stderr, "failure in lu decomposition\n");
					exit(1);
				}
			}
		}

		/* parameter */
		for(i=0; i<param_num; i++){
			if(param[i]->depending_rule != NULL && param[i]->depending_rule->is_assignment){
				param[i]->temp_value = param[i]->k[step];
			}else{
				for(l=0; l<step; l++){
					param[i]->temp_value = param[i]->value + param[i]->k[l]*dt*rk_ce[step-1][l];
				}
				error = calc_by_algebraic(algEq, cycle, dt, *reverse_time, time, res, print_interval, err_zero_flag);
				if (error == 1) {
					fprintf(stderr, "failure in lu decomposition\n");
					exit(1);
				}
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
						for(l=0; l<step; l++){
							comp[i]->including_species[j]->temp_value = comp[i]->including_species[j]->temp_value*comp[i]->temp_value/(comp[i]->value + comp[i]->k[l]*dt*rk_ce[step-1][l]);
						}
					}
				}
				/* new code end */
				for(l=0; l<step; l++){
					comp[i]->temp_value = comp[i]->value + comp[i]->k[l]*dt*rk_ce[step-1][l];
				}
				error = calc_by_algebraic(algEq, cycle, dt, *reverse_time, time, res, print_interval, err_zero_flag);
			}
		}
		/* species reference */
		for(i=0; i<spr_num; i++){
			if(spr[i]->depending_rule != NULL && spr[i]->depending_rule->is_assignment){
				spr[i]->temp_value = spr[i]->k[step];
			}else{
				for(l=0; l<step; l++){
					spr[i]->temp_value = spr[i]->value + spr[i]->k[l]*dt*rk_ce[step-1][l];
				}
				error = calc_by_algebraic(algEq, cycle, dt, *reverse_time, time, res, print_interval, err_zero_flag);
				if (error == 1) {
					fprintf(stderr, "failure in lu decomposition\n");
					exit(1);
				}
			}
		}
    }
  }
}


/* calculate the temp_value after the stage is changed */
int calc_by_algebraic(myAlgebraicEquations *algEq, int cycle, double dt, double  reverse_time, double* time, myResult* result, int print_interval, int* err_zero_flag){
	int error;
	unsigned int i,j;
	double* constant_vector = NULL;
	double **coefficient_matrix = NULL;
	int* alg_pivot = NULL;

    /* calc temp value algebraic by algebraic */
    if(algEq != NULL){
		if(algEq->num_of_algebraic_variables > 1){
			/* initialize pivot */
			for(i=0; i<algEq->num_of_algebraic_variables; i++){
				alg_pivot[i] = i;
			}
			for(i=0; i<algEq->num_of_algebraic_variables; i++){
				for(j=0; j<algEq->num_of_algebraic_variables; j++){
					coefficient_matrix[i][j] = calcf(algEq->coefficient_matrix[i][j], dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
				}
			}
			for(i=0; i<algEq->num_of_algebraic_variables; i++){
				constant_vector[i] = -calcf(algEq->constant_vector[i], dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
			}
			/* LU decompostion */
			error = lu_decomposition(coefficient_matrix, alg_pivot, algEq->num_of_algebraic_variables);
			if(error == 0){/* failure in LU decomposition */

				return 1;
			}
			/* forward substitution & backward substitution */
			lu_solve(coefficient_matrix, alg_pivot, algEq->num_of_algebraic_variables, constant_vector);
			for(i=0; i<algEq->num_of_alg_target_sp; i++){
				algEq->alg_target_species[i]->target_species->temp_value = constant_vector[algEq->alg_target_species[i]->order];
			}
			for(i=0; i<algEq->num_of_alg_target_param; i++){
				algEq->alg_target_parameter[i]->target_parameter->temp_value = constant_vector[algEq->alg_target_parameter[i]->order];
			}
			for(i=0; i<algEq->num_of_alg_target_comp; i++){
				/* new code */
				for(j=0; j<algEq->alg_target_compartment[i]->target_compartment->num_of_including_species; j++){
					if(algEq->alg_target_compartment[i]->target_compartment->including_species[j]->is_concentration){
						algEq->alg_target_compartment[i]->target_compartment->including_species[j]->temp_value = algEq->alg_target_compartment[i]->target_compartment->including_species[j]->temp_value*algEq->alg_target_compartment[i]->target_compartment->temp_value/constant_vector[algEq->alg_target_compartment[i]->order];
					}
				}
				/* new code end  */
				algEq->alg_target_compartment[i]->target_compartment->temp_value = constant_vector[algEq->alg_target_compartment[i]->order];
			}
		}else{
			if(algEq->target_species != NULL){
				algEq->target_species->temp_value = -calcf(algEq->constant, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag)/calcf(algEq->coefficient, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
			}
			if(algEq->target_parameter != NULL){
				algEq->target_parameter->temp_value = -calcf(algEq->constant, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag)/calcf(algEq->coefficient, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
			}
			if(algEq->target_compartment != NULL){
				/* new code */
			for(i=0; i<algEq->target_compartment->num_of_including_species; i++){
				if(algEq->target_compartment->including_species[i]->is_concentration){
					algEq->target_compartment->including_species[i]->temp_value = algEq->target_compartment->including_species[i]->temp_value*algEq->target_compartment->temp_value/(-calcf(algEq->constant, dt, cycle, &reverse_time, 0, time, time,  result, print_interval, err_zero_flag)/calcf(algEq->coefficient, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag));
				}
			}
			/* new code end */
			algEq->target_compartment->temp_value = -calcf(algEq->constant, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag)/calcf(algEq->coefficient, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
			}
		}
    }
	return 0;
}

/* for variable stepsize, calculate the temp_value after the stage is changed */

void calc_by_assignment(mySpecies *sp[], unsigned int sp_num, myParameter *param[], unsigned int param_num, myCompartment *comp[], unsigned int comp_num, mySpeciesReference *spr[], unsigned int spr_num, double dt, int cycle, double reverse_time, double* time, myResult* result, int print_interval, int* err_zero_flag){
	unsigned int i;
	/* calc temp value by assignment */
	for(i=0; i<sp_num; i++){
		if(sp[i]->depending_rule != NULL && sp[i]->depending_rule->is_assignment){
			sp[i]->temp_value = calcf(sp[i]->depending_rule->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
		}
	}
	for(i=0; i<param_num; i++){
		if(param[i]->depending_rule != NULL && param[i]->depending_rule->is_assignment){
			param[i]->temp_value = calcf(param[i]->depending_rule->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
		}
	}
	for(i=0; i<comp_num; i++){
		if(comp[i]->depending_rule != NULL && comp[i]->depending_rule->is_assignment){
			comp[i]->temp_value = calcf(comp[i]->depending_rule->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
		}
	}
	for(i=0; i<spr_num; i++){
		if(spr[i]->depending_rule != NULL && spr[i]->depending_rule->is_assignment){
			spr[i]->temp_value = calcf(spr[i]->depending_rule->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
		}
	}
}
