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
#include <time.h>
#include "../libsbmlsim/libsbmlsim.h"

/* coefficient matrix for Adams Bashforth */
double c_e[4][4] = {{1.0, 0, 0, 0}, /* AB1 (Euler) : order = 0 */
  {3.0/2.0, -1.0/2.0, 0, 0}, /* AB2 : order = 1 */
  {23.0/12.0, -16.0/12.0, 5.0/12.0, 0}, /* AB3 : order = 2 */
  {55.0/24.0, -59.0/24.0, 37.0/24.0, -9.0/24.0}}; /* AB4 : order = 3 */

double calc_explicit_formula(int order, double k1, double k2, double k3, double k4){
  return c_e[order][0]*k1 + c_e[order][1]*k2 + c_e[order][2]*k3 + c_e[order][3]*k4;
}

myResult* simulate_explicit(Model_t *m, myResult* result, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], myRule *rule[], myEvent *event[], myInitialAssignment *initAssign[], myAlgebraicEquations *algEq, timeVariantAssignments *timeVarAssign, double sim_time, double dt, int print_interval, double *time, int order, int print_amount, allocated_memory *mem){
  unsigned int i, j;
  int cycle;
  int error;
  int end_cycle = get_end_cycle(sim_time, dt);
  double reverse_time;
  double *value_time_p = result->values_time;
  double *value_sp_p = result->values_sp;
  double *value_param_p = result->values_param;
  double *value_comp_p = result->values_comp;
  double **coefficient_matrix = NULL;
  double *constant_vector = NULL;
  int *alg_pivot = NULL;
  double reactants_numerator, products_numerator;
  double min_value;
  double *init_val;

  /* num of SBase objects */
  unsigned int num_of_species = Model_getNumSpecies(m);
  unsigned int num_of_parameters = Model_getNumParameters(m);
  unsigned int num_of_compartments = Model_getNumCompartments(m);
  unsigned int num_of_reactions = Model_getNumReactions(m);
  unsigned int num_of_rules = Model_getNumRules(m);
  unsigned int num_of_events = Model_getNumEvents(m);
  unsigned int num_of_initialAssignments = Model_getNumInitialAssignments(m);
  /* num of variables whose quantity is not a constant */
  unsigned int num_of_all_var_species = 0;
  unsigned int num_of_all_var_parameters = 0;
  unsigned int num_of_all_var_compartments = 0;
  unsigned int num_of_all_var_species_reference = 0;
  /* num of variables (which is NOT changed by assignment nor algebraic rule) */
  unsigned int num_of_var_species = 0;
  unsigned int num_of_var_parameters = 0;
  unsigned int num_of_var_compartments = 0;
  unsigned int num_of_var_species_reference = 0;
  /* All variables (whose quantity is not a constant) */
  mySpecies **all_var_sp;           /* all variable species */
  myParameter **all_var_param;      /* all variable parameters */
  myCompartment **all_var_comp;     /* all variable compartments */
  mySpeciesReference **all_var_spr; /* all varialbe SpeciesReferences */
  /* variables (which is NOT changed by assignment nor algebraic rule) */
  mySpecies **var_sp; 
  myParameter **var_param;
  myCompartment **var_comp;
  mySpeciesReference **var_spr;

  set_seed();

  check_num(num_of_species, num_of_parameters, num_of_compartments, num_of_reactions, &num_of_all_var_species, &num_of_all_var_parameters, &num_of_all_var_compartments, &num_of_all_var_species_reference, &num_of_var_species, &num_of_var_parameters, &num_of_var_compartments, &num_of_var_species_reference, sp, param, comp, re);

  /* create objects */
  all_var_sp = (mySpecies **)malloc(sizeof(mySpecies *) * num_of_all_var_species);
  all_var_param = (myParameter **)malloc(sizeof(myParameter *) * num_of_all_var_parameters);
  all_var_comp = (myCompartment **)malloc(sizeof(myCompartment *) * num_of_all_var_compartments);
  all_var_spr = (mySpeciesReference **)malloc(sizeof(mySpeciesReference *) * num_of_all_var_species_reference);
  var_sp = (mySpecies **)malloc(sizeof(mySpecies *) * num_of_var_species);
  var_param = (myParameter **)malloc(sizeof(myParameter *) * num_of_var_parameters);
  var_comp = (myCompartment **)malloc(sizeof(myCompartment *) * num_of_var_compartments);
  var_spr = (mySpeciesReference **)malloc(sizeof(mySpeciesReference *) * num_of_var_species_reference);
  /* mySpecies *all_var_sp[num_of_all_var_species]; */
  /* myParameter *all_var_param[num_of_all_var_parameters]; */
  /* myCompartment *all_var_comp[num_of_all_var_compartments]; */
  /* mySpeciesReference *all_var_spr[num_of_all_var_species_reference]; */
  /* mySpecies *var_sp[num_of_var_species]; */
  /* myParameter *var_param[num_of_var_parameters]; */
  /* myCompartment *var_comp[num_of_var_compartments]; */
  /* mySpeciesReference *var_spr[num_of_var_species_reference]; */

  create_calc_object_list(num_of_species, num_of_parameters, num_of_compartments, num_of_reactions, all_var_sp, all_var_param, all_var_comp, all_var_spr, var_sp, var_param, var_comp, var_spr, sp, param, comp, re);

  if(algEq != NULL){
    coefficient_matrix = (double**)malloc(sizeof(double*)*(algEq->num_of_algebraic_variables));
    for(i=0; i<algEq->num_of_algebraic_variables; i++){
      coefficient_matrix[i] = (double*)malloc(sizeof(double)*(algEq->num_of_algebraic_variables));
    }
    constant_vector = (double*)malloc(sizeof(double)*(algEq->num_of_algebraic_variables));
    alg_pivot = (int*)malloc(sizeof(int)*(algEq->num_of_algebraic_variables));
  }

  PRG_TRACE(("Simulation for [%s] Starts!\n", Model_getId(m)));
  cycle = 0;

  /* initialize delay_val */
  initialize_delay_val(sp, num_of_species, param, num_of_parameters, comp, num_of_compartments, re, num_of_reactions, sim_time, dt, 0);

  /* calc temp value by assignment */
  for(i=0; i<num_of_all_var_species; i++){
    if(all_var_sp[i]->depending_rule != NULL && all_var_sp[i]->depending_rule->is_assignment){
      all_var_sp[i]->temp_value = calc(all_var_sp[i]->depending_rule->eq, dt, cycle, &reverse_time, 0);
    }
  }
  for(i=0; i<num_of_all_var_parameters; i++){
    if(all_var_param[i]->depending_rule != NULL && all_var_param[i]->depending_rule->is_assignment){
      all_var_param[i]->temp_value = calc(all_var_param[i]->depending_rule->eq, dt, cycle, &reverse_time, 0);
    }
  }
  for(i=0; i<num_of_all_var_compartments; i++){
    if(all_var_comp[i]->depending_rule != NULL && all_var_comp[i]->depending_rule->is_assignment){
      all_var_comp[i]->temp_value = calc(all_var_comp[i]->depending_rule->eq, dt, cycle, &reverse_time, 0);
    }
  }
  for(i=0; i<num_of_all_var_species_reference; i++){
    if(all_var_spr[i]->depending_rule != NULL && all_var_spr[i]->depending_rule->is_assignment){
      all_var_spr[i]->temp_value = calc(all_var_spr[i]->depending_rule->eq, dt, cycle, &reverse_time, 0);
    }
  }
  /* forwarding value */
  forwarding_value(all_var_sp, num_of_all_var_species, all_var_param, num_of_all_var_parameters, all_var_comp, num_of_all_var_compartments, all_var_spr, num_of_all_var_species_reference);

  /* initialize delay_val */
  initialize_delay_val(sp, num_of_species, param, num_of_parameters, comp, num_of_compartments, re, num_of_reactions, sim_time, dt, 0);

  /* calc InitialAssignment */
  calc_initial_assignment(initAssign, num_of_initialAssignments, dt, cycle, &reverse_time);

  /* initialize delay_val */
  initialize_delay_val(sp, num_of_species, param, num_of_parameters, comp, num_of_compartments, re, num_of_reactions, sim_time, dt, 0);

  /* rewriting for explicit delay */
  for(i=0; i<num_of_initialAssignments; i++){
    for(j=0; j<initAssign[i]->eq->math_length; j++){
      if(initAssign[i]->eq->number[j] == time){
        TRACE(("time is replaced with reverse time\n"));
        initAssign[i]->eq->number[j] = &reverse_time;
      }else if(initAssign[i]->eq->number[j] != NULL){
        init_val = (double*)malloc(sizeof(double));
        *init_val = *initAssign[i]->eq->number[j];
        mem->memory[mem->num_of_allocated_memory++] = init_val;
        initAssign[i]->eq->number[j] = init_val;
      }
    }
  }
  for(i=0; i<timeVarAssign->num_of_time_variant_assignments; i++){
    for(j=0; j<timeVarAssign->eq[i]->math_length; j++){
      if(timeVarAssign->eq[i]->number[j] == time){
        TRACE(("time is replaced with reverse time\n"));
        timeVarAssign->eq[i]->number[j] = &reverse_time;
      }else if(timeVarAssign->eq[i]->number[j] != NULL){
        init_val = (double*)malloc(sizeof(double));
        *init_val = *timeVarAssign->eq[i]->number[j];
        mem->memory[mem->num_of_allocated_memory++] = init_val;
        timeVarAssign->eq[i]->number[j] = init_val;
      }
    }
  }

  /* initialize delay_val */
  initialize_delay_val(sp, num_of_species, param, num_of_parameters, comp, num_of_compartments, re, num_of_reactions, sim_time, dt, 0);

  /* calc temp value by assignment */
  for(i=0; i<num_of_all_var_species; i++){
    if(all_var_sp[i]->depending_rule != NULL && all_var_sp[i]->depending_rule->is_assignment){
      all_var_sp[i]->temp_value = calc(all_var_sp[i]->depending_rule->eq, dt, cycle, &reverse_time, 0);
    }
  }
  for(i=0; i<num_of_all_var_parameters; i++){
    if(all_var_param[i]->depending_rule != NULL && all_var_param[i]->depending_rule->is_assignment){
      all_var_param[i]->temp_value = calc(all_var_param[i]->depending_rule->eq, dt, cycle, &reverse_time, 0);
    }
  }
  for(i=0; i<num_of_all_var_compartments; i++){
    if(all_var_comp[i]->depending_rule != NULL && all_var_comp[i]->depending_rule->is_assignment){
      all_var_comp[i]->temp_value = calc(all_var_comp[i]->depending_rule->eq, dt, cycle, &reverse_time, 0);
    }
  }
  for(i=0; i<num_of_all_var_species_reference; i++){
    if(all_var_spr[i]->depending_rule != NULL && all_var_spr[i]->depending_rule->is_assignment){
      all_var_spr[i]->temp_value = calc(all_var_spr[i]->depending_rule->eq, dt, cycle, &reverse_time, 0);
    }
  }
  /* forwarding value */
  forwarding_value(all_var_sp, num_of_all_var_species, all_var_param, num_of_all_var_parameters, all_var_comp, num_of_all_var_compartments, all_var_spr, num_of_all_var_species_reference);

  /* initialize delay_val */
  initialize_delay_val(sp, num_of_species, param, num_of_parameters, comp, num_of_compartments, re, num_of_reactions, sim_time, dt, 0);

  /* calc temp value algebraic by algebraic */
  if(algEq != NULL){
    if(algEq->num_of_algebraic_variables > 1){
      /* initialize pivot */
      for(i=0; i<algEq->num_of_algebraic_variables; i++){
        alg_pivot[i] = i;
      }
      for(i=0; i<algEq->num_of_algebraic_variables; i++){
        for(j=0; j<algEq->num_of_algebraic_variables; j++){
          coefficient_matrix[i][j] = calc(algEq->coefficient_matrix[i][j], dt, cycle, &reverse_time, 0);
          /* TRACE(("coefficient matrix[%d][%d] = %lf\n", i, j, coefficient_matrix[i][j])); */
        }
      }
      for(i=0; i<algEq->num_of_algebraic_variables; i++){
        constant_vector[i] = -calc(algEq->constant_vector[i], dt, cycle, &reverse_time, 0);
        /* TRACE(("constant vector[%d] = %lf\n", i, constant_vector[i])); */
      }
      /* LU decompostion */
      error = lu_decomposition(coefficient_matrix, alg_pivot, algEq->num_of_algebraic_variables);
      if(error == 0){/* failure in LU decomposition */
        return NULL;
      }
      /* forward substitution & backward substitution */
      lu_solve(coefficient_matrix, alg_pivot, algEq->num_of_algebraic_variables, constant_vector);
      /*       for(i=0; i<algEq->num_of_algebraic_variables; i++){ */
      /* 	TRACE(("ans[%d] = %lf\n", i, constant_vector[i])); */
      /*       } */
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
       /* new code end */
        algEq->alg_target_compartment[i]->target_compartment->temp_value = constant_vector[algEq->alg_target_compartment[i]->order];
      }    
    }else{
      if(algEq->target_species != NULL){
        algEq->target_species->temp_value = -calc(algEq->constant, dt, cycle, &reverse_time, 0)/calc(algEq->coefficient, dt, cycle, &reverse_time, 0);
      }
      if(algEq->target_parameter != NULL){
        algEq->target_parameter->temp_value = -calc(algEq->constant, dt, cycle, &reverse_time, 0)/calc(algEq->coefficient, dt, cycle, &reverse_time, 0);
      }
      if(algEq->target_compartment != NULL){
        /* new code */
        for(i=0; i<algEq->target_compartment->num_of_including_species; i++){
          if(algEq->target_compartment->including_species[i]->is_concentration){
            algEq->target_compartment->including_species[i]->temp_value = algEq->target_compartment->including_species[i]->temp_value*algEq->target_compartment->temp_value/(-calc(algEq->constant, dt, cycle, &reverse_time, 0)/calc(algEq->coefficient, dt, cycle, &reverse_time, 0));
          }
        }
       /* new code end */
        algEq->target_compartment->temp_value = -calc(algEq->constant, dt, cycle, &reverse_time, 0)/calc(algEq->coefficient, dt, cycle, &reverse_time, 0);
      }
    }
    /* forwarding value */
    forwarding_value(all_var_sp, num_of_all_var_species, all_var_param, num_of_all_var_parameters, all_var_comp, num_of_all_var_compartments, all_var_spr, num_of_all_var_species_reference);
  }

  /* initialize delay_val */
  initialize_delay_val(sp, num_of_species, param, num_of_parameters, comp, num_of_compartments, re, num_of_reactions, sim_time, dt, 1);

  /* cycle start */
  for(cycle=0; cycle<=end_cycle; cycle++){
    /* calculate unreversible fast reaction */
    for(i=0; i<num_of_reactions; i++){
      if(re[i]->is_fast && !re[i]->is_reversible){
        if(calc(re[i]->eq, dt, cycle, &reverse_time, 0) > 0){
          min_value = DBL_MAX;
          for(j=0; j<re[i]->num_of_reactants; j++){
            if(min_value > re[i]->reactants[j]->mySp->value/calc(re[i]->reactants[j]->eq, dt, cycle, &reverse_time, 0)){
              min_value = re[i]->reactants[j]->mySp->value/calc(re[i]->reactants[j]->eq, dt, cycle, &reverse_time, 0);
            }
          }
          for(j=0; j<re[i]->num_of_products; j++){
            if(!Species_getBoundaryCondition(re[i]->products[j]->mySp->origin)){
              re[i]->products[j]->mySp->value += calc(re[i]->products[j]->eq, dt, cycle, &reverse_time, 0)*min_value;
              re[i]->products[j]->mySp->temp_value = re[i]->products[j]->mySp->value;
            }
          }
          for(j=0; j<re[i]->num_of_reactants; j++){
            if(!Species_getBoundaryCondition(re[i]->reactants[j]->mySp->origin)){
              re[i]->reactants[j]->mySp->value -= calc(re[i]->reactants[j]->eq, dt, cycle, &reverse_time, 0)*min_value;
              re[i]->reactants[j]->mySp->temp_value = re[i]->reactants[j]->mySp->value;
            }
          }
        }
      }
    }
    /* calculate reversible fast reaction */
    for(i=0; i<num_of_reactions; i++){
      if(re[i]->is_fast && re[i]->is_reversible){
        if(!(Species_getBoundaryCondition(re[i]->products[0]->mySp->origin) 
              && Species_getBoundaryCondition(re[i]->reactants[0]->mySp->origin))){
          products_numerator = calc(re[i]->products_equili_numerator, dt, cycle, &reverse_time, 0);
          reactants_numerator = calc(re[i]->reactants_equili_numerator, dt, cycle, &reverse_time, 0);
          if(products_numerator > 0 || reactants_numerator > 0){
            if(Species_getBoundaryCondition(re[i]->products[0]->mySp->origin)){
              re[i]->reactants[0]->mySp->value = (reactants_numerator/products_numerator)*re[i]->products[0]->mySp->value;
              re[i]->reactants[0]->mySp->temp_value = re[i]->reactants[0]->mySp->value;
            }else if(Species_getBoundaryCondition(re[i]->reactants[0]->mySp->origin)){
              re[i]->products[0]->mySp->value = (products_numerator/reactants_numerator)*re[i]->reactants[0]->mySp->value;
              re[i]->products[0]->mySp->temp_value = re[i]->products[0]->mySp->value;	    
            }else{
              re[i]->products[0]->mySp->value = (products_numerator/(products_numerator+reactants_numerator))*(re[i]->products[0]->mySp->temp_value+re[i]->reactants[0]->mySp->temp_value);
              re[i]->reactants[0]->mySp->value = (reactants_numerator/(products_numerator+reactants_numerator))*(re[i]->products[0]->mySp->temp_value+re[i]->reactants[0]->mySp->temp_value);
              re[i]->products[0]->mySp->temp_value = re[i]->products[0]->mySp->value;
              re[i]->reactants[0]->mySp->temp_value = re[i]->reactants[0]->mySp->value;
            }
          }
        }
      }
    }

    /* event */
    calc_event(event, num_of_events, dt, *time, cycle, &reverse_time);

    /* substitute delay val */
    substitute_delay_val(sp, num_of_species, param, num_of_parameters, comp, num_of_compartments, re, num_of_reactions, cycle);

    /* progress */
    if(cycle%(int)(end_cycle/10) == 0){
      PRG_TRACE(("%3d %%\n", (int)(100*((double)cycle/(double)end_cycle))));
      PRG_TRACE(("\x1b[1A"));
      PRG_TRACE(("\x1b[5D"));
    }
    /* print result */
	if(cycle%print_interval == 0) {
      /*  Time */
      *value_time_p = *time;
      value_time_p++;
      /*  Species */
      for(i=0; i<num_of_species; i++){
        if(print_amount){
          if(sp[i]->is_concentration){
            *value_sp_p = sp[i]->value*sp[i]->locating_compartment->value;
          }else{
            *value_sp_p = sp[i]->value;
          }
        }else{
          if(sp[i]->is_amount){
            *value_sp_p = sp[i]->value/sp[i]->locating_compartment->value;
          }else{
            *value_sp_p = sp[i]->value;
          }
        }
        value_sp_p++;
      }
      /*  Parameter */
      for(i=0; i<num_of_parameters; i++){
        *value_param_p = param[i]->value;
        value_param_p++;
      }
      /*  Compartment */
      for(i=0; i<num_of_compartments; i++){
        *value_comp_p = comp[i]->value;
        value_comp_p++;
      }
    }

    /* time increase */
    *time = (cycle+1)*dt;

    if(order == 4){/* runge kutta       */
      calc_k(all_var_sp, num_of_all_var_species, all_var_param, num_of_all_var_parameters, all_var_comp, num_of_all_var_compartments, all_var_spr, num_of_all_var_species_reference, re, num_of_reactions, rule, num_of_rules, cycle, dt, &reverse_time, 1, 1);
      calc_temp_value(all_var_sp, num_of_all_var_species, all_var_param, num_of_all_var_parameters, all_var_comp, num_of_all_var_compartments, all_var_spr, num_of_all_var_species_reference, dt, 1);      
    } else {/* Adams-Bashforth */
      /* calc k */
      calc_k(var_sp, num_of_var_species, var_param, num_of_var_parameters, var_comp, num_of_var_compartments, var_spr, num_of_var_species_reference, re, num_of_reactions, rule, num_of_rules, cycle, dt, &reverse_time, 0, 1);      
      /* calc temp value by Adams Bashforth */
      for(i=0; i<num_of_var_species; i++){
        var_sp[i]->temp_value = var_sp[i]->value + calc_explicit_formula(order, var_sp[i]->k[0], var_sp[i]->prev_k[0], var_sp[i]->prev_k[1], var_sp[i]->prev_k[2])*dt;
      }
      for(i=0; i<num_of_var_parameters; i++){
        var_param[i]->temp_value = var_param[i]->value + calc_explicit_formula(order, var_param[i]->k[0], var_param[i]->prev_k[0], var_param[i]->prev_k[1], var_param[i]->prev_k[2])*dt;
      }
      for(i=0; i<num_of_var_compartments; i++){
        var_comp[i]->temp_value = var_comp[i]->value + calc_explicit_formula(order, var_comp[i]->k[0], var_comp[i]->prev_k[0], var_comp[i]->prev_k[1], var_comp[i]->prev_k[2])*dt;
      }
      for(i=0; i<num_of_var_species_reference; i++){
        var_spr[i]->temp_value = var_spr[i]->value + calc_explicit_formula(order, var_spr[i]->k[0], var_spr[i]->prev_k[0], var_spr[i]->prev_k[1], var_spr[i]->prev_k[2])*dt;
      }      
      /* calc temp value by assignment */
      for(i=0; i<num_of_all_var_species; i++){
        if(all_var_sp[i]->depending_rule != NULL && all_var_sp[i]->depending_rule->is_assignment){
          all_var_sp[i]->temp_value = calc(all_var_sp[i]->depending_rule->eq, dt, cycle, &reverse_time, 0);
        }
      }
      for(i=0; i<num_of_all_var_parameters; i++){
        if(all_var_param[i]->depending_rule != NULL && all_var_param[i]->depending_rule->is_assignment){
          all_var_param[i]->temp_value = calc(all_var_param[i]->depending_rule->eq, dt, cycle, &reverse_time, 0);
        }
      }
      for(i=0; i<num_of_all_var_compartments; i++){
        if(all_var_comp[i]->depending_rule != NULL && all_var_comp[i]->depending_rule->is_assignment){
          all_var_comp[i]->temp_value = calc(all_var_comp[i]->depending_rule->eq, dt, cycle, &reverse_time, 0);
        }
      }
      for(i=0; i<num_of_all_var_species_reference; i++){
        if(all_var_spr[i]->depending_rule != NULL && all_var_spr[i]->depending_rule->is_assignment){
          all_var_spr[i]->temp_value = calc(all_var_spr[i]->depending_rule->eq, dt, cycle, &reverse_time, 0);
        }
      }
    }

    /* calc temp value algebraic by algebraic */
    if(algEq != NULL){
      if(algEq->num_of_algebraic_variables > 1){
        /* initialize pivot */
        for(i=0; i<algEq->num_of_algebraic_variables; i++){
          alg_pivot[i] = i;
        }
        for(i=0; i<algEq->num_of_algebraic_variables; i++){
          for(j=0; j<algEq->num_of_algebraic_variables; j++){
            coefficient_matrix[i][j] = calc(algEq->coefficient_matrix[i][j], dt, cycle, &reverse_time, 0);
          }
        }
        for(i=0; i<algEq->num_of_algebraic_variables; i++){
          constant_vector[i] = -calc(algEq->constant_vector[i], dt, cycle, &reverse_time, 0);
        }
        /* LU decompostion */
        error = lu_decomposition(coefficient_matrix, alg_pivot, algEq->num_of_algebraic_variables);
        if(error == 0){/* failure in LU decomposition */
          return NULL;
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
          algEq->target_species->temp_value = -calc(algEq->constant, dt, cycle, &reverse_time, 0)/calc(algEq->coefficient, dt, cycle, &reverse_time, 0);
        }
        if(algEq->target_parameter != NULL){
          algEq->target_parameter->temp_value = -calc(algEq->constant, dt, cycle, &reverse_time, 0)/calc(algEq->coefficient, dt, cycle, &reverse_time, 0);
        }
        if(algEq->target_compartment != NULL){
          /* new code */
          for(i=0; i<algEq->target_compartment->num_of_including_species; i++){
            if(algEq->target_compartment->including_species[i]->is_concentration){
              algEq->target_compartment->including_species[i]->temp_value = algEq->target_compartment->including_species[i]->temp_value*algEq->target_compartment->temp_value/(-calc(algEq->constant, dt, cycle, &reverse_time, 0)/calc(algEq->coefficient, dt, cycle, &reverse_time, 0));
            }
          }
         /* new code end */
          algEq->target_compartment->temp_value = -calc(algEq->constant, dt, cycle, &reverse_time, 0)/calc(algEq->coefficient, dt, cycle, &reverse_time, 0);
        }
      }
    }

    /* preserve prev_value and prev_k */
    for(i=0; i<num_of_var_species; i++){
      var_sp[i]->prev_val[2] = var_sp[i]->prev_val[1];
      var_sp[i]->prev_val[1] = var_sp[i]->prev_val[0];
      var_sp[i]->prev_val[0] = var_sp[i]->value;
      var_sp[i]->prev_k[2] = var_sp[i]->prev_k[1];
      var_sp[i]->prev_k[1] = var_sp[i]->prev_k[0];
      var_sp[i]->prev_k[0] = var_sp[i]->k[0];
    }
    for(i=0; i<num_of_var_parameters; i++){
      var_param[i]->prev_val[2] = var_param[i]->prev_val[1];
      var_param[i]->prev_val[1] = var_param[i]->prev_val[0];
      var_param[i]->prev_val[0] = var_param[i]->value;
      var_param[i]->prev_k[2] = var_param[i]->prev_k[1];
      var_param[i]->prev_k[1] = var_param[i]->prev_k[0];
      var_param[i]->prev_k[0] = var_param[i]->k[0];
    }
    for(i=0; i<num_of_var_compartments; i++){
      var_comp[i]->prev_val[2] = var_comp[i]->prev_val[1];
      var_comp[i]->prev_val[1] = var_comp[i]->prev_val[0];
      var_comp[i]->prev_val[0] = var_comp[i]->value;
      var_comp[i]->prev_k[2] = var_comp[i]->prev_k[1];
      var_comp[i]->prev_k[1] = var_comp[i]->prev_k[0];
      var_comp[i]->prev_k[0] = var_comp[i]->k[0];
    }
    for(i=0; i<num_of_var_species_reference; i++){
      var_spr[i]->prev_val[2] = var_spr[i]->prev_val[1];
      var_spr[i]->prev_val[1] = var_spr[i]->prev_val[0];
      var_spr[i]->prev_val[0] = var_spr[i]->value;
      var_spr[i]->prev_k[2] = var_spr[i]->prev_k[1];
      var_spr[i]->prev_k[1] = var_spr[i]->prev_k[0];
      var_spr[i]->prev_k[0] = var_spr[i]->k[0];
    }
    /* forwarding value */
    forwarding_value(all_var_sp, num_of_all_var_species, all_var_param, num_of_all_var_parameters, all_var_comp, num_of_all_var_compartments, all_var_spr, num_of_all_var_species_reference);
    /* time increase */
	//*time += dt;
  }
  PRG_TRACE(("Simulation for [%s] Ends!\n", Model_getId(m)));
  if(algEq != NULL){
	  for(i=0; i<algEq->num_of_algebraic_variables; i++){
		  free(coefficient_matrix[i]);
	  }
	  free(coefficient_matrix);
	  free(constant_vector);
	  free(alg_pivot);
  }
  free(all_var_sp);
  free(all_var_param);
  free(all_var_comp);
  free(all_var_spr);
  free(var_sp);
  free(var_param);
  free(var_comp);
  free(var_spr);
  return result;
}

myResult* simulate_explicitf(Model_t *m, myResult* result, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], myRule *rule[], myEvent *event[], myInitialAssignment *initAssign[], myAlgebraicEquations *algEq, timeVariantAssignments *timeVarAssign, double sim_time, double dt, int print_interval, double *time, int order, int print_amount, allocated_memory *mem, double atol, double rtol, double facmax, copied_AST *cp_AST, int* err_zero_flag){
  unsigned int i, j;
  int cycle;
  int error;
  int end_cycle = get_end_cycle(sim_time, dt);
  double reverse_time;
  double *value_time_p = result->values_time;
  int model_has_delay = 0;
  double* value_time_p_fordelay;

  double *value_sp_p = result->values_sp;
  double *value_param_p = result->values_param;
  double *value_comp_p = result->values_comp;
  double **coefficient_matrix = NULL;
  double *constant_vector = NULL;
  int *alg_pivot = NULL;
  double reactants_numerator, products_numerator;
  double min_value;
  double *init_val;

  /*for variable step-size integration*/
  double sum_error = 0.0;
  int max_index = end_cycle + 1;
  int max_result_index = end_cycle * print_interval + 1;
  int err_zero_flag2 = 0;
  int ode_num = 0;
  double fixed_dt = dt;
  double* value_fixed_time_p = (double*)malloc(sizeof(double) * (end_cycle + 2));
  double* value_fixed_time_h1 = value_fixed_time_p;
  double* value_fixed_time_h2 = value_fixed_time_p;
  int* ode_check;

  /* num of SBase objects */
  unsigned int num_of_species = Model_getNumSpecies(m);
  unsigned int num_of_parameters = Model_getNumParameters(m);
  unsigned int num_of_compartments = Model_getNumCompartments(m);
  unsigned int num_of_reactions = Model_getNumReactions(m);
  unsigned int num_of_rules = Model_getNumRules(m);
  unsigned int num_of_events = Model_getNumEvents(m);
  unsigned int num_of_initialAssignments = Model_getNumInitialAssignments(m);
  /* num of variables whose quantity is not a constant */
  unsigned int num_of_all_var_species = 0;
  unsigned int num_of_all_var_parameters = 0;
  unsigned int num_of_all_var_compartments = 0;
  unsigned int num_of_all_var_species_reference = 0;
  /* num of variables (which is NOT changed by assignment nor algebraic rule) */
  unsigned int num_of_var_species = 0;
  unsigned int num_of_var_parameters = 0;
  unsigned int num_of_var_compartments = 0;
  unsigned int num_of_var_species_reference = 0;
  /* All variables (whose quantity is not a constant) */
  mySpecies **all_var_sp;           /* all variable species */
  myParameter **all_var_param;      /* all variable parameters */
  myCompartment **all_var_comp;     /* all variable compartments */
  mySpeciesReference **all_var_spr; /* all varialbe SpeciesReferences */
  /* variables (which is NOT changed by assignment nor algebraic rule) */
  mySpecies **var_sp;
  myParameter **var_param;
  myCompartment **var_comp;
  mySpeciesReference **var_spr;
  /* adjustment of tolerance*/
  /* atol *= 1e-10; */
  /* rtol *= 1e-07; */

  *(err_zero_flag) = 0;
  /* make fixed time step array */
  for(cycle=0; cycle<=end_cycle+1; cycle++){
	  *value_fixed_time_p = cycle * dt;
	  value_fixed_time_p++;
  }

  /* to count the number of ODE */
  ode_check = (int*)malloc(sizeof(int) * num_of_species);
  for(i=0; i<num_of_species; i++){
	  *(ode_check + i) = 0;
  }

  set_seed();

  check_num(num_of_species, num_of_parameters, num_of_compartments, num_of_reactions, &num_of_all_var_species, &num_of_all_var_parameters, &num_of_all_var_compartments, &num_of_all_var_species_reference, &num_of_var_species, &num_of_var_parameters, &num_of_var_compartments, &num_of_var_species_reference, sp, param, comp, re);

  /* create objects */
  all_var_sp = (mySpecies **)malloc(sizeof(mySpecies *) * num_of_all_var_species);
  all_var_param = (myParameter **)malloc(sizeof(myParameter *) * num_of_all_var_parameters);
  all_var_comp = (myCompartment **)malloc(sizeof(myCompartment *) * num_of_all_var_compartments);
  all_var_spr = (mySpeciesReference **)malloc(sizeof(mySpeciesReference *) * num_of_all_var_species_reference);
  var_sp = (mySpecies **)malloc(sizeof(mySpecies *) * num_of_var_species);
  var_param = (myParameter **)malloc(sizeof(myParameter *) * num_of_var_parameters);
  var_comp = (myCompartment **)malloc(sizeof(myCompartment *) * num_of_var_compartments);
  var_spr = (mySpeciesReference **)malloc(sizeof(mySpeciesReference *) * num_of_var_species_reference);
  /* mySpecies *all_var_sp[num_of_all_var_species]; */
  /* myParameter *all_var_param[num_of_all_var_parameters]; */
  /* myCompartment *all_var_comp[num_of_all_var_compartments]; */
  /* mySpeciesReference *all_var_spr[num_of_all_var_species_reference]; */
  /* mySpecies *var_sp[num_of_var_species]; */
  /* myParameter *var_param[num_of_var_parameters]; */
  /* myCompartment *var_comp[num_of_var_compartments]; */
  /* mySpeciesReference *var_spr[num_of_var_species_reference]; */

  create_calc_object_list(num_of_species, num_of_parameters, num_of_compartments, num_of_reactions, all_var_sp, all_var_param, all_var_comp, all_var_spr, var_sp, var_param, var_comp, var_spr, sp, param, comp, re);

  if(algEq != NULL){
    coefficient_matrix = (double**)malloc(sizeof(double*)*(algEq->num_of_algebraic_variables));
    for(i=0; i<algEq->num_of_algebraic_variables; i++){
      coefficient_matrix[i] = (double*)malloc(sizeof(double)*(algEq->num_of_algebraic_variables));
    }
    constant_vector = (double*)malloc(sizeof(double)*(algEq->num_of_algebraic_variables));
    alg_pivot = (int*)malloc(sizeof(int)*(algEq->num_of_algebraic_variables));
  }

  PRG_TRACE(("Simulation for [%s] Starts!\n", Model_getId(m)));
  cycle = 0;
  *(time) = 0.0;

  /* initialize delay_val */
  initialize_delay_valf(sp, num_of_species, param, num_of_parameters, comp, num_of_compartments, re, num_of_reactions, sim_time, dt, 0);

  /* new code */
  /* if model has delay, save the value of time at all step
	 (for linear approximation) */
  for(i=0; i<num_of_species; i++){
	  if(sp[i]->delay_val != NULL){
		  model_has_delay = 1;
	  }
  }
  for(i=0; i<num_of_parameters; i++){
	  if(param[i]->delay_val != NULL){
		  model_has_delay = 1;
	  }
  }
  for(i=0; i<num_of_compartments; i++){
	  if(comp[i]->delay_val != NULL){
		  model_has_delay = 1;
	  }
  }
  if (model_has_delay){
	  result->num_of_delay_rows = result->num_of_rows;
	  result->values_time_fordelay = (double*)malloc(sizeof(double) * max_result_index);
	  value_time_p_fordelay = result->values_time_fordelay;
  }
  /* new code end */

  /* calc temp value by assignment */
  for(i=0; i<num_of_all_var_species; i++){
    if(all_var_sp[i]->depending_rule != NULL && all_var_sp[i]->depending_rule->is_assignment){
		all_var_sp[i]->temp_value = calcf(all_var_sp[i]->depending_rule->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
    }
  }
  for(i=0; i<num_of_all_var_parameters; i++){
    if(all_var_param[i]->depending_rule != NULL && all_var_param[i]->depending_rule->is_assignment){
		all_var_param[i]->temp_value = calcf(all_var_param[i]->depending_rule->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
    }
  }
  for(i=0; i<num_of_all_var_compartments; i++){
    if(all_var_comp[i]->depending_rule != NULL && all_var_comp[i]->depending_rule->is_assignment){
		all_var_comp[i]->temp_value = calcf(all_var_comp[i]->depending_rule->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
    }
  }
  for(i=0; i<num_of_all_var_species_reference; i++){
    if(all_var_spr[i]->depending_rule != NULL && all_var_spr[i]->depending_rule->is_assignment){
		all_var_spr[i]->temp_value = calcf(all_var_spr[i]->depending_rule->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
    }
  }
    /* forwarding value */
    forwarding_value(all_var_sp, num_of_all_var_species, all_var_param, num_of_all_var_parameters, all_var_comp, num_of_all_var_compartments, all_var_spr, num_of_all_var_species_reference);

  /* initialize delay_val */
  initialize_delay_valf(sp, num_of_species, param, num_of_parameters, comp, num_of_compartments, re, num_of_reactions, sim_time, dt, 0);

  /* calc InitialAssignment */
  calc_initial_assignmentf(initAssign, num_of_initialAssignments, dt, cycle, &reverse_time, time, result, print_interval, err_zero_flag);

  /* initialize delay_val */
  initialize_delay_valf(sp, num_of_species, param, num_of_parameters, comp, num_of_compartments, re, num_of_reactions, sim_time, dt, 0);

  /* rewriting for explicit delay */
  for(i=0; i<num_of_initialAssignments; i++){
    for(j=0; j<initAssign[i]->eq->math_length; j++){
		if(initAssign[i]->eq->op[j] == AST_NAME_TIME){
			TRACE(("time is replaced with reverse time\n"));
			initAssign[i]->eq->time_reverse_flag = 1;
			initAssign[i]->eq->reverse_time = &reverse_time;
		}else if(initAssign[i]->eq->number[j] != NULL){
			init_val = (double*)malloc(sizeof(double));
			*init_val = *initAssign[i]->eq->number[j];
			mem->memory[mem->num_of_allocated_memory++] = init_val;
			initAssign[i]->eq->number[j] = init_val;
  }
    }
  }
  for(i=0; i<timeVarAssign->num_of_time_variant_assignments; i++){
	  for(j=0; j<timeVarAssign->eq[i]->math_length; j++){
		  if(timeVarAssign->eq[i]->op[j] == AST_NAME_TIME){
			  TRACE(("time is replaced with reverse time\n"));
			  timeVarAssign->eq[i]->time_reverse_flag = 1;
			  timeVarAssign->eq[i]->reverse_time = &reverse_time;
		  }else if(timeVarAssign->eq[i]->number[j] != NULL){
			  init_val = (double*)malloc(sizeof(double));
			  *init_val = *timeVarAssign->eq[i]->number[j];
			  mem->memory[mem->num_of_allocated_memory++] = init_val;
			  timeVarAssign->eq[i]->number[j] = init_val;
		  }
	  }
  }
  /* initialize delay_val */
  initialize_delay_valf(sp, num_of_species, param, num_of_parameters, comp, num_of_compartments, re, num_of_reactions, sim_time, dt, 0);

  /* calc temp value by assignment */
  for(i=0; i<num_of_all_var_species; i++){
	  if(all_var_sp[i]->depending_rule != NULL && all_var_sp[i]->depending_rule->is_assignment){
		  all_var_sp[i]->temp_value = calcf(all_var_sp[i]->depending_rule->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
	  }
  }
  for(i=0; i<num_of_all_var_parameters; i++){
	  if(all_var_param[i]->depending_rule != NULL && all_var_param[i]->depending_rule->is_assignment){
		  all_var_param[i]->temp_value = calcf(all_var_param[i]->depending_rule->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
	  }
  }
  for(i=0; i<num_of_all_var_compartments; i++){
	  if(all_var_comp[i]->depending_rule != NULL && all_var_comp[i]->depending_rule->is_assignment){
		  all_var_comp[i]->temp_value = calcf(all_var_comp[i]->depending_rule->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
	  }
  }
  for(i=0; i<num_of_all_var_species_reference; i++){
	  if(all_var_spr[i]->depending_rule != NULL && all_var_spr[i]->depending_rule->is_assignment){
		  all_var_spr[i]->temp_value = calcf(all_var_spr[i]->depending_rule->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
	  }
  }

  /* forwarding value */
  forwarding_value(all_var_sp, num_of_all_var_species, all_var_param, num_of_all_var_parameters, all_var_comp, num_of_all_var_compartments, all_var_spr, num_of_all_var_species_reference);

  /* initialize delay_val */
  initialize_delay_valf(sp, num_of_species, param, num_of_parameters, comp, num_of_compartments, re, num_of_reactions, sim_time, dt, 0);

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
          /* TRACE(("coefficient matrix[%d][%d] = %lf\n", i, j, coefficient_matrix[i][j])); */
        }
      }
      for(i=0; i<algEq->num_of_algebraic_variables; i++){
		  constant_vector[i] = -calcf(algEq->constant_vector[i], dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
        /* TRACE(("constant vector[%d] = %lf\n", i, constant_vector[i])); */
      }
      /* LU decompostion */
      error = lu_decomposition(coefficient_matrix, alg_pivot, algEq->num_of_algebraic_variables);
      if(error == 0){/* failure in LU decomposition */
        return NULL;
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
		  /* new code end */
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
					algEq->target_compartment->including_species[i]->temp_value = algEq->target_compartment->including_species[i]->temp_value*algEq->target_compartment->temp_value/(-calcf(algEq->constant, dt, cycle, &reverse_time, 0,  time, time, result, print_interval, err_zero_flag)/calcf(algEq->coefficient, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag));
				}
			}
			/* new code end */
			algEq->target_compartment->temp_value = -calcf(algEq->constant, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag)/calcf(algEq->coefficient, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
		}
    }
    /* forwarding value */
    forwarding_value(all_var_sp, num_of_all_var_species, all_var_param, num_of_all_var_parameters, all_var_comp, num_of_all_var_compartments, all_var_spr, num_of_all_var_species_reference);
  }

  /* initialize delay_val */
  initialize_delay_valf(sp, num_of_species, param, num_of_parameters, comp, num_of_compartments, re, num_of_reactions, sim_time, dt, 1);

  /* for variable stepsize, calculate number of ODE */
  /* rate rules of Species */
  for(i=0; i<num_of_all_var_species; i++){
	  if(all_var_sp[i]->depending_rule != NULL && all_var_sp[i]->depending_rule->is_rate){
		  ode_num++;
	  }
  }
  /* rate rules of parameters */
  for(i=0; i<num_of_all_var_parameters; i++){
	  if(all_var_param[i]->depending_rule != NULL && all_var_param[i]->depending_rule->is_rate){
		  ode_num++;
	  }
  }
  /* rate rules of compartments */
  for(i=0; i<num_of_all_var_compartments; i++){
	  if(all_var_comp[i]->depending_rule != NULL && all_var_comp[i]->depending_rule->is_rate){
		  ode_num++;
	  }
  }
  /* reactions */
  for(i=0; i<num_of_reactions; i++){
      if(!re[i]->is_fast){
		  for(j=0; j<re[i]->num_of_products; j++){
			  if(!Species_getBoundaryCondition(re[i]->products[j]->mySp->origin)){
				  if(count_ode(sp, num_of_species, ode_check, re[i]->products[j]->mySp->origin)){
				  ode_num++;
				  }
			  }
		  }
		  for(j=0; j<re[i]->num_of_reactants; j++){
			  if(!Species_getBoundaryCondition(re[i]->reactants[j]->mySp->origin)){
				  if(count_ode(sp, num_of_species, ode_check, re[i]->reactants[j]->mySp->origin)){
				  ode_num++;
				  }
			  }
		  }
      }
  }

  /*  --- cycle start --- */
  while(1) {
	  /* if necessary, reallocate objects */
	  if (cycle >= max_index && *(err_zero_flag) == 0) {
		  realloc_mySBML_objects(m, sp, num_of_species, param, num_of_parameters, comp, num_of_compartments, re, num_of_reactions, rule, event, num_of_events, initAssign, &timeVarAssign, cp_AST, sim_time, max_index);
		  max_index = (max_index - 1) * 2 + 1;
	  }
	  /* if necessary, reallocate result contents */
	  if (cycle >= max_result_index && *(err_zero_flag) == 0 && model_has_delay) {
		  reallocate_result_objects(result, &value_time_p_fordelay, max_result_index);
		  max_result_index *= 2;
	  }

	  /* calculate unreversible fast reaction */
	  for(i=0; i<num_of_reactions; i++){
		if(re[i]->is_fast && !re[i]->is_reversible){
			if(calcf(re[i]->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag) > 0){
				min_value = DBL_MAX;
				for(j=0; j<re[i]->num_of_reactants; j++){
					if(min_value > re[i]->reactants[j]->mySp->value/calcf(re[i]->reactants[j]->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag)){
						min_value = re[i]->reactants[j]->mySp->value/calcf(re[i]->reactants[j]->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
					}
				}
				for(j=0; j<re[i]->num_of_products; j++){
					if(!Species_getBoundaryCondition(re[i]->products[j]->mySp->origin)){
						re[i]->products[j]->mySp->value += calcf(re[i]->products[j]->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag)*min_value;
						re[i]->products[j]->mySp->temp_value = re[i]->products[j]->mySp->value;
					}
				}
				for(j=0; j<re[i]->num_of_reactants; j++){
					if(!Species_getBoundaryCondition(re[i]->reactants[j]->mySp->origin)){
						re[i]->reactants[j]->mySp->value -= calcf(re[i]->reactants[j]->eq, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag)*min_value;
						re[i]->reactants[j]->mySp->temp_value = re[i]->reactants[j]->mySp->value;
					}
				}
			}
		}
	  }

	  /* calculate reversible fast reaction */
	  for(i=0; i<num_of_reactions; i++){
		  if(re[i]->is_fast && re[i]->is_reversible){
			  if(!(Species_getBoundaryCondition(re[i]->products[0]->mySp->origin)
				   && Species_getBoundaryCondition(re[i]->reactants[0]->mySp->origin))){
				  products_numerator = calcf(re[i]->products_equili_numerator, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
				  reactants_numerator = calcf(re[i]->reactants_equili_numerator, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
				  if(products_numerator > 0 || reactants_numerator > 0){
					  if(Species_getBoundaryCondition(re[i]->products[0]->mySp->origin)){
						  re[i]->reactants[0]->mySp->value = (reactants_numerator/products_numerator)*re[i]->products[0]->mySp->value;
						  re[i]->reactants[0]->mySp->temp_value = re[i]->reactants[0]->mySp->value;
					  }else if(Species_getBoundaryCondition(re[i]->reactants[0]->mySp->origin)){
						  re[i]->products[0]->mySp->value = (products_numerator/reactants_numerator)*re[i]->reactants[0]->mySp->value;
						  re[i]->products[0]->mySp->temp_value = re[i]->products[0]->mySp->value;
					  }else{
						  re[i]->products[0]->mySp->value = (products_numerator/(products_numerator+reactants_numerator))*(re[i]->products[0]->mySp->temp_value+re[i]->reactants[0]->mySp->temp_value);
						  re[i]->reactants[0]->mySp->value = (reactants_numerator/(products_numerator+reactants_numerator))*(re[i]->products[0]->mySp->temp_value+re[i]->reactants[0]->mySp->temp_value);
						  re[i]->products[0]->mySp->temp_value = re[i]->products[0]->mySp->value;
						  re[i]->reactants[0]->mySp->temp_value = re[i]->reactants[0]->mySp->value;
					}
				  }
			  }
		  }
	  }

	  /* event */
	  calc_eventf(event, num_of_events, dt, *time, cycle, &reverse_time, result, print_interval, err_zero_flag);

	  /* substitute delay val */
	  substitute_delay_valf(sp, num_of_species, param, num_of_parameters, comp, num_of_compartments, re, num_of_reactions, cycle);

	  /* print result for fixed stepsize
		 (only when cycle is multiple of print_interval)　*/
	  if (*(err_zero_flag) == 1) {
		  if (cycle % print_interval == 0) {
			  /*  Time */
			  *value_time_p = *time;
			  value_time_p++;
			  /*  Species */
			  for(i=0; i<num_of_species; i++){
				  if(print_amount){
					  if(sp[i]->is_concentration){
						  *value_sp_p = sp[i]->value*sp[i]->locating_compartment->value;
					  }else{
						  *value_sp_p = sp[i]->value;
					  }
				  }else{
					  if(sp[i]->is_amount){
						  *value_sp_p = sp[i]->value/sp[i]->locating_compartment->value;
					  }else{
					  *value_sp_p = sp[i]->value;
					  }
				  }
				  value_sp_p++;
			  }
			  /*  Parameter */
			  for(i=0; i<num_of_parameters; i++){
				  *value_param_p = param[i]->value;
				  value_param_p++;
			  }
			  /*  Compartment */
			  for(i=0; i<num_of_compartments; i++){
				  *value_comp_p = comp[i]->value;
				  value_comp_p++;
			  }
			  /* new code */
			  /* save the time values for delay */
			  if(model_has_delay){
				  *value_time_p_fordelay = *time;
				  value_time_p_fordelay++;
			  }
			  /* new code end */
		  }
	  }

	  /* print result for variable stepsize*/
	  if (*(err_zero_flag) == 0 && cycle == 0){
		  /*  Time */
		  *value_time_p = *value_fixed_time_h1;
		  value_time_p++;
		  value_fixed_time_h1++;
		  /*  Species */
		  for(i=0; i<num_of_species; i++){
			  if(print_amount){
				  if(sp[i]->is_concentration){
					  *value_sp_p = sp[i]->value*sp[i]->locating_compartment->value;
				  }else{
					  *value_sp_p = sp[i]->value;
				  }
			  }else{
				  if(sp[i]->is_amount){
					  *value_sp_p = sp[i]->value/sp[i]->locating_compartment->value;
				  }else{
					  *value_sp_p = sp[i]->value;
				  }
			  }
			  value_sp_p++;
		  }
		  /*  Parameter */
		  for(i=0; i<num_of_parameters; i++){
			  *value_param_p = param[i]->value;
			  value_param_p++;
		  }
		  /*  Compartment */
		  for(i=0; i<num_of_compartments; i++){
			  *value_comp_p = comp[i]->value;
			  value_comp_p++;
		  }
		  /* new code */
		  /* save the all time values for delay */
		  if(model_has_delay){
			  *value_time_p_fordelay = *time;
			  value_time_p_fordelay++;
		  }
		  /* new code end */
	  }

	  if (order == 5 || order == 6){/* Runge-Kutta-Fehlberg or Cash-Karp */
		  do {
			  calc_kf(all_var_sp, num_of_all_var_species, all_var_param, num_of_all_var_parameters, all_var_comp, num_of_all_var_compartments, all_var_spr, num_of_all_var_species_reference, re, num_of_reactions, rule, num_of_rules, cycle, dt, &reverse_time, 1, 1, time, result, algEq, print_interval, err_zero_flag, order);
			  sum_error = calc_sum_error(all_var_sp, num_of_all_var_species, all_var_param, num_of_all_var_parameters, all_var_comp, num_of_all_var_compartments, all_var_spr, num_of_all_var_species_reference, dt, 1, atol, rtol, &ode_num, 0, order);
			  if((sum_error == 0.0 && ode_num == 0 )||
				 (sum_error == 0.0 && cycle == 0) ||
				 *(err_zero_flag) == 1) {
				  if (*(err_zero_flag) == 0) {
					  *(err_zero_flag) = 1;
					  dt = fixed_dt / print_interval;
					  if (sum_error == 0.0 && ode_num == 0){
						  TRACE(("no ode model \n"));
					  }else{
						  TRACE(("no error model \n"));
					  }
				  }
				  break;
			  }else if (sum_error == 0.0 || err_zero_flag2){
				  if (!err_zero_flag2) {
					  TRACE(("in some intervals of integration, the model has no error\n"));
					  err_zero_flag2 = 1;
					  dt = fixed_dt / print_interval;
				  }
				  break;
			  }else{
				  dt *= my_fmin(facmax, my_fmax(0.9 * pow((1.0/sum_error), 1.0/(order+1)), 1e-04));/* order should be 5 */
			  }
			  if (*(time) == *(time) + dt){
				  dt = calc_eps(dt);
				  break;
			  }
        /*
			  if (!isfinite(dt) || !isfinite(sum_error)) {
				  TRACE(("stepsize ERROR\n"));
				  //write_csv(result, "variable.csv"); save the result data by variable step-size
				  exit(1);
          }
          */
		  } while (sum_error > 1);
		  /* time increase */
		  if (*(err_zero_flag) == 0) {
			  *time += dt;
		  }else {
			  *(time) = (cycle + 1) * dt;
		  }
		  sum_error = calc_sum_error(all_var_sp, num_of_all_var_species, all_var_param, num_of_all_var_parameters, all_var_comp, num_of_all_var_compartments, all_var_spr, num_of_all_var_species_reference, dt, 1, atol, rtol, &ode_num, 1, order);
		  /* new code */
		  calc_by_assignment(all_var_sp, num_of_all_var_species, all_var_param, num_of_all_var_parameters, all_var_comp, num_of_all_var_compartments, all_var_spr, num_of_all_var_species_reference, dt, cycle, reverse_time, time, result, print_interval, err_zero_flag);
		  /* new code end */
	  }

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
				return NULL;
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
					algEq->target_compartment->including_species[i]->temp_value = algEq->target_compartment->including_species[i]->temp_value*algEq->target_compartment->temp_value/(-calcf(algEq->constant, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag)/calcf(algEq->coefficient, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag));
				}
			}
			/* new code end */
			algEq->target_compartment->temp_value = -calcf(algEq->constant, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag)/calcf(algEq->coefficient, dt, cycle, &reverse_time, 0, time, time, result, print_interval, err_zero_flag);
			}
		}
    }

    /* preserve prev_value and prev_k (for multistep solution)*/
    for(i=0; i<num_of_var_species; i++){
		var_sp[i]->prev_val[2] = var_sp[i]->prev_val[1];
		var_sp[i]->prev_val[1] = var_sp[i]->prev_val[0];
		var_sp[i]->prev_val[0] = var_sp[i]->value;
		var_sp[i]->prev_k[2] = var_sp[i]->prev_k[1];
		var_sp[i]->prev_k[1] = var_sp[i]->prev_k[0];
		var_sp[i]->prev_k[0] = var_sp[i]->k[0];
    }
    for(i=0; i<num_of_var_parameters; i++){
		var_param[i]->prev_val[2] = var_param[i]->prev_val[1];
		var_param[i]->prev_val[1] = var_param[i]->prev_val[0];
		var_param[i]->prev_val[0] = var_param[i]->value;
		var_param[i]->prev_k[2] = var_param[i]->prev_k[1];
		var_param[i]->prev_k[1] = var_param[i]->prev_k[0];
		var_param[i]->prev_k[0] = var_param[i]->k[0];
    }
    for(i=0; i<num_of_var_compartments; i++){
		var_comp[i]->prev_val[2] = var_comp[i]->prev_val[1];
		var_comp[i]->prev_val[1] = var_comp[i]->prev_val[0];
		var_comp[i]->prev_val[0] = var_comp[i]->value;
		var_comp[i]->prev_k[2] = var_comp[i]->prev_k[1];
		var_comp[i]->prev_k[1] = var_comp[i]->prev_k[0];
		var_comp[i]->prev_k[0] = var_comp[i]->k[0];
    }
    for(i=0; i<num_of_var_species_reference; i++){
		var_spr[i]->prev_val[2] = var_spr[i]->prev_val[1];
		var_spr[i]->prev_val[1] = var_spr[i]->prev_val[0];
		var_spr[i]->prev_val[0] = var_spr[i]->value;
		var_spr[i]->prev_k[2] = var_spr[i]->prev_k[1];
		var_spr[i]->prev_k[1] = var_spr[i]->prev_k[0];
		var_spr[i]->prev_k[0] = var_spr[i]->k[0];
    }

	/* print result for variable stepsize*/
	if (*(err_zero_flag) == 0){
		/* new code */
		/* save the all time values for delay */
		if(model_has_delay){
			*value_time_p_fordelay = *time;
			value_time_p_fordelay++;
		}
		/* new code end */

		if(*time >= *value_fixed_time_h1){
			while(1){
				/*  Species */
				for(i=0; i<num_of_species; i++){
					if(print_amount){
						if(sp[i]->is_concentration){
							*value_sp_p = approximate_printresult_linearly(sp[i]->value, sp[i]->temp_value, *time, *time - dt, *value_fixed_time_h1)*approximate_printresult_linearly(sp[i]->locating_compartment->value, sp[i]->locating_compartment->temp_value, *time, *time - dt, *value_fixed_time_h1);
						}else{
							*value_sp_p = approximate_printresult_linearly(sp[i]->value, sp[i]->temp_value, *time, *time - dt, *value_fixed_time_h1);
						}
					}else{
						if(sp[i]->is_amount){
							*value_sp_p = approximate_printresult_linearly(sp[i]->value, sp[i]->temp_value, *time, *time - dt, *value_fixed_time_h1)/approximate_printresult_linearly(sp[i]->locating_compartment->value, sp[i]->locating_compartment->temp_value, *time, *time - dt, *value_fixed_time_h1);
						}else{
							*value_sp_p = approximate_printresult_linearly(sp[i]->value, sp[i]->temp_value, *time, *time - dt, *value_fixed_time_h1);
						}
					}
					value_sp_p++;
				}
				  /*  Parameter */
				for(i=0; i<num_of_parameters; i++){
					*value_param_p = approximate_printresult_linearly(param[i]->value, param[i]->temp_value, *time, *time - dt, *value_fixed_time_h1);
					value_param_p++;
				}
				/*  Compartment */
				for(i=0; i<num_of_compartments; i++){
					*value_comp_p = approximate_printresult_linearly(comp[i]->value, comp[i]->temp_value, *time, *time - dt, *value_fixed_time_h1);
					value_comp_p++;
				}
				/*  Time */
				*value_time_p = *value_fixed_time_h1;
				value_time_p++;
				value_fixed_time_h1++;

				if (*value_fixed_time_h1 > sim_time || *time < *value_fixed_time_h1){
					break;
				}
			}
		}
	}
    /* forwarding value */
    forwarding_value(all_var_sp, num_of_all_var_species, all_var_param, num_of_all_var_parameters, all_var_comp, num_of_all_var_compartments, all_var_spr, num_of_all_var_species_reference);

    /* cycle increase */
	cycle++;
	/* --- branching for simulation end --- */
	/* for variable stepsize */
	if (*(err_zero_flag) == 0 && *time >= sim_time) {
		break;
	}
	/* for fixed stepsize */
	if (*(err_zero_flag) == 1 && cycle > end_cycle * print_interval) {
		break;
	}
  }

  PRG_TRACE(("Simulation for [%s] Ends!\n", Model_getId(m)));

  if (model_has_delay){
    free(result->values_time_fordelay);
  }
  free(ode_check);
  free(value_fixed_time_h2);
  if(algEq != NULL){
    for(i=0; i<algEq->num_of_algebraic_variables; i++){
      free(coefficient_matrix[i]);
    }
    free(coefficient_matrix);
    free(constant_vector);
    free(alg_pivot);
  }
  free(all_var_sp);
  free(all_var_param);
  free(all_var_comp);
  free(all_var_spr);
  free(var_sp);
  free(var_param);
  free(var_comp);
  free(var_spr);
  //print_result(result);
  return result;
}

