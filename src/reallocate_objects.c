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

void reallocate_objects(Model_t *m, mySpecies *sp[], unsigned int num_of_species, myParameter *param[], unsigned int num_of_parameters, myCompartment *comp[], unsigned int num_of_compartments, myReaction *re[], unsigned int num_of_reactions, myRule *rule[], myEvent *ev[], unsigned int num_of_events, myInitialAssignment *myInitAssign[], timeVariantAssignments **timeVarAssign, copied_AST *cp_AST, double sim_time, int max_index) {
	unsigned int i, j, l;
	int k;
	int new_max_index = (max_index - 1) * 2 + 1;
	/* for substitution of reallocated address */
	int flag = 0;
	ASTNode_t* node;
	ASTNode_t *times_node, *conv_factor_node;
	unsigned int math_length = 0;
	/* 1) reallocate delay value arrays */
	/* species */
	for(i=0; i<num_of_species; i++){
		if(sp[i]->delay_val != NULL){
			sp[i]->delay_val = (double**)realloc(sp[i]->delay_val, sizeof(double*) * new_max_index);
			for(k=0; k<=new_max_index; k++) {
				if (k <max_index) {
					sp[i]->delay_val[k] = (double *)realloc(sp[i]->delay_val[k], sizeof(double) * 6);
				}else {
					sp[i]->delay_val[k] = (double *)malloc(sizeof(double) * 6);
				}
			}
			flag = 1;
		}
	}
	/* parameter */
	for(i=0; i<num_of_parameters; i++){
		if(param[i]->delay_val != NULL){
			param[i]->delay_val = (double**)realloc(param[i]->delay_val, sizeof(double*) * new_max_index);
			for(k=0; k<=new_max_index; k++) {
				if (k <max_index) {
					param[i]->delay_val[k] = (double *)realloc(param[i]->delay_val[k], sizeof(double) * 6);
				}else {
					param[i]->delay_val[k] = (double *)malloc(sizeof(double) * 6);
				}
			}
			flag = 1;
		}
	}
	/* compartment*/
	for(i=0; i<num_of_compartments; i++){
		if(comp[i]->delay_val != NULL){
			comp[i]->delay_val = (double**)realloc(comp[i]->delay_val, sizeof(double*) * new_max_index );
			for(k=0; k<=new_max_index; k++) {
				if (k <max_index) {
					comp[i]->delay_val[k] = (double *)realloc(comp[i]->delay_val[k], sizeof(double) * 6);
				}else {
					comp[i]->delay_val[k] = (double *)malloc(sizeof(double) * 6);
				}
			}
			flag = 1;
		}
	}
	/* products and reactants */
	for(i=0; i<num_of_reactions; i++){
		for(j=0; j<re[i]->num_of_products; j++){
			if(re[i]->products[j]->delay_val != NULL){
				re[i]->products[j]->delay_val = (double**)realloc(re[i]->products[j]->delay_val, sizeof(double*)* new_max_index);
				for(k=0; k<=new_max_index; k++) {
					if (k <max_index) {
						re[i]->products[j]->delay_val[k] = (double *)realloc(re[i]->products[j]->delay_val[k], sizeof(double) * 6);
					}else {
						re[i]->products[j]->delay_val[k] = (double *)malloc(sizeof(double) * 6);
					}
				}
				flag = 1;
			}
		}
		for(j=0; j<re[i]->num_of_reactants; j++){
			if(re[i]->reactants[j]->delay_val != NULL){
				re[i]->reactants[j]->delay_val = (double**)realloc(re[i]->reactants[j]->delay_val, sizeof(double*)* new_max_index);
				for(k=0; k<=new_max_index; k++) {
					if (k <max_index) {
						re[i]->reactants[j]->delay_val[k] = (double *)realloc(re[i]->reactants[j]->delay_val[k], sizeof(double) * 6);
					}else {
						re[i]->reactants[j]->delay_val[k] = (double *)malloc(sizeof(double) * 6);
					}
				}
				flag = 1;
			}
		}
	}
	if(flag) {
		/* connect delay_val with all the equation->delay_number */
		/* initial assignments */
		for(i=0; i<Model_getNumInitialAssignments(m); i++){
			node = (ASTNode_t*)InitialAssignment_getMath(myInitAssign[i]->origin);
			node = ASTNode_deepCopy(node);
			check_AST(node, NULL);
			/* alter_tree_structure */
			alter_tree_structure(m, &node, NULL, 0, cp_AST);
			math_length = connect_delayval_with_eq(m, myInitAssign[i]->eq, sp, param, comp, re, node, 0);
		}
		/* time var assignments */
		(*timeVarAssign)->num_of_time_variant_assignments = 0;
		for(i=0; i<Model_getNumRules(m); i++){
			node = (ASTNode_t*)Rule_getMath(rule[i]->origin);
			if(Rule_isAssignment(rule[i]->origin) && include_time(node, 0)){
				node = ASTNode_deepCopy(node);
				check_AST(node, NULL);
				alter_tree_structure(m, &node, NULL, 0, cp_AST);
				for(j=0; j<Model_getNumSpecies(m); j++){
					if(strcmp(Species_getId(sp[j]->origin), (char*)Rule_getVariable(rule[i]->origin)) == 0){
						if(sp[j]->is_amount
						   && !sp[j]->has_only_substance_units
						   && Compartment_getSpatialDimensions(sp[j]->locating_compartment->origin) != 0){
							assignment_alter_tree_structure(&node, (char*)Compartment_getId(sp[j]->locating_compartment->origin), 0);
						}else if(sp[j]->is_concentration
								 && sp[j]->has_only_substance_units){
							assignment_alter_tree_structure(&node, (char*)Compartment_getId(sp[j]->locating_compartment->origin), 1);
						}
						break;
					}
				}
				check_AST(node, NULL);
				math_length = connect_delayval_with_eq(m, (*timeVarAssign)->eq[(*timeVarAssign)->num_of_time_variant_assignments], sp, param, comp, re, node, 0);
			}
		}
		/* reaction equation */
		for(i=0; i<num_of_reactions; i++){
			node = (ASTNode_t*)KineticLaw_getMath(Reaction_getKineticLaw(re[i]->origin));
			node = ASTNode_deepCopy(node);
			check_AST(node, NULL);
			alter_tree_structure(m, &node, NULL, 0, cp_AST);
			set_local_para_as_value(node, Reaction_getKineticLaw(re[i]->origin));
			check_AST(node, NULL);
			math_length = connect_delayval_with_eq(m, re[i]->eq, sp, param, comp, re, node, 0);
			/* products */
			for(j=0; j<re[i]->num_of_products; j++){
				node = (ASTNode_t*)StoichiometryMath_getMath(SpeciesReference_getStoichiometryMath(re[i]->products[j]->origin));
				if(node != NULL){/* l2v4 */
					node = ASTNode_deepCopy(node);
				}else if(SpeciesReference_isSetStoichiometry(re[i]->products[j]->origin) && !SpeciesReference_isSetId(re[i]->products[j]->origin)){
					node = ASTNode_createWithType(AST_REAL);
					ASTNode_setReal(node, SpeciesReference_getStoichiometry(re[i]->products[j]->origin));
				}else if(SpeciesReference_isSetStoichiometry(re[i]->products[j]->origin) && SpeciesReference_isSetId(re[i]->products[j]->origin)){
					node = ASTNode_createWithType(AST_NAME);
					ASTNode_setName(node, SpeciesReference_getId(re[i]->products[j]->origin));
				}
				if(node == NULL){
					if(my_isnan(SpeciesReference_getStoichiometry(re[i]->products[j]->origin)) && SpeciesReference_isSetId(re[i]->products[j]->origin)){/* l3v1 */
						node = ASTNode_createWithType(AST_NAME);
						ASTNode_setName(node, SpeciesReference_getId(re[i]->products[j]->origin));
					}else if(my_isnan(SpeciesReference_getStoichiometry(re[i]->products[j]->origin))){/* l3v1 */
						node = ASTNode_createWithType(AST_REAL);
						ASTNode_setReal(node, 1.0);
					}else{/* l2v4 */
						node = ASTNode_createWithType(AST_REAL);
						ASTNode_setReal(node, SpeciesReference_getStoichiometry(re[i]->products[j]->origin));
					}
				}
				check_AST(node, NULL);
				alter_tree_structure(m, &node, NULL, 0, cp_AST);
				/* unit */
				if(re[i]->products[j]->mySp->is_concentration){
					assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(re[i]->products[j]->mySp->origin))), 1);
				}
				/* unit */
				/* conversion factor */
				if(Species_isSetConversionFactor(re[i]->products[j]->mySp->origin)){
					times_node = ASTNode_createWithType(AST_TIMES);
					conv_factor_node = ASTNode_createWithType(AST_NAME);
					ASTNode_setName(conv_factor_node, Species_getConversionFactor(re[i]->products[j]->mySp->origin));
					ASTNode_addChild(times_node, conv_factor_node);
					ASTNode_addChild(times_node, node);
					node = times_node;
				}else if(Model_isSetConversionFactor(m)){
					times_node = ASTNode_createWithType(AST_TIMES);
					conv_factor_node = ASTNode_createWithType(AST_NAME);
					ASTNode_setName(conv_factor_node, Model_getConversionFactor(m));
					ASTNode_addChild(times_node, conv_factor_node);
					ASTNode_addChild(times_node, node);
					node = times_node;
				}
				/* conversion factor */
				math_length = connect_delayval_with_eq(m, re[i]->products[j]->eq, sp, param, comp, re, node, 0);
				check_AST(node, NULL);
			}
			/* reactants */
			for(j=0; j<re[i]->num_of_reactants; j++){
				node = (ASTNode_t*)StoichiometryMath_getMath(SpeciesReference_getStoichiometryMath(re[i]->reactants[j]->origin));
				if(node != NULL){/* l2v4 */
					node = ASTNode_deepCopy(node);
				}else if(SpeciesReference_isSetStoichiometry(re[i]->reactants[j]->origin) && !SpeciesReference_isSetId(re[i]->reactants[j]->origin)){
					node = ASTNode_createWithType(AST_REAL);
					ASTNode_setReal(node, SpeciesReference_getStoichiometry(re[i]->reactants[j]->origin));
				}else if(SpeciesReference_isSetStoichiometry(re[i]->reactants[j]->origin) && SpeciesReference_isSetId(re[i]->reactants[j]->origin)){
					node = ASTNode_createWithType(AST_NAME);
					ASTNode_setName(node, SpeciesReference_getId(re[i]->reactants[j]->origin));
				}
				if(node == NULL){
					if(my_isnan(SpeciesReference_getStoichiometry(re[i]->reactants[j]->origin)) && SpeciesReference_isSetId(re[i]->reactants[j]->origin)){/* l3v1 */
						node = ASTNode_createWithType(AST_NAME);
						ASTNode_setName(node, SpeciesReference_getId(re[i]->reactants[j]->origin));
					}else if(my_isnan(SpeciesReference_getStoichiometry(re[i]->reactants[j]->origin))){/* l3v1 */
						node = ASTNode_createWithType(AST_REAL);
						ASTNode_setReal(node, 1.0);
					}else{/* l2v4 */
						node = ASTNode_createWithType(AST_REAL);
						ASTNode_setReal(node, SpeciesReference_getStoichiometry(re[i]->reactants[j]->origin));
					}
				}
				check_AST(node, NULL);
				alter_tree_structure(m, &node, NULL, 0, cp_AST);
				/* unit */
				if(re[i]->reactants[j]->mySp->is_concentration){
					assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(re[i]->reactants[j]->mySp->origin))), 1);
				}
				/* unit */
				/* conversion factor */
				if(Species_isSetConversionFactor(re[i]->reactants[j]->mySp->origin)){
					times_node = ASTNode_createWithType(AST_TIMES);
					conv_factor_node = ASTNode_createWithType(AST_NAME);
					ASTNode_setName(conv_factor_node, Species_getConversionFactor(re[i]->reactants[j]->mySp->origin));
					ASTNode_addChild(times_node, conv_factor_node);
					ASTNode_addChild(times_node, node);
					node = times_node;
				}else if(Model_isSetConversionFactor(m)){
					times_node = ASTNode_createWithType(AST_TIMES);
					conv_factor_node = ASTNode_createWithType(AST_NAME);
					ASTNode_setName(conv_factor_node, Model_getConversionFactor(m));
					ASTNode_addChild(times_node, conv_factor_node);
					ASTNode_addChild(times_node, node);
					node = times_node;
				}
				/* conversion factor */
				math_length =connect_delayval_with_eq(m, re[i]->reactants[j]->eq, sp, param, comp, re, node, 0);
				check_AST(node, NULL);
			}
		}

		/* create myRules */
		for(i=0; i<Model_getNumRules(m); i++){
			if(Rule_isRate(rule[i]->origin) || Rule_isAssignment(rule[i]->origin)){
				flag = 1;
				for(j=0; j<num_of_species; j++){
					if(strcmp(Rule_getVariable(rule[i]->origin), Species_getId(sp[j]->origin)) == 0){
						flag = 0;
						break;
					}
				}
				if(flag){
					for(j=0; j<num_of_parameters; j++){
						if(strcmp(Rule_getVariable(rule[i]->origin), Parameter_getId(param[j]->origin)) == 0){
							flag = 0;
							break;
						}
					}
				}
				if(flag){
					for(j=0; j<num_of_compartments; j++){
						if(strcmp(Rule_getVariable(rule[i]->origin), Compartment_getId(comp[j]->origin)) == 0){
							flag = 0;
							break;
						}
					}
				}
				if(flag){
					for(j=0; j<num_of_reactions; j++){
						for(l=0; l<re[j]->num_of_products; l++){
							if(SpeciesReference_isSetId(re[j]->products[l]->origin)
							   && strcmp(Rule_getVariable(rule[i]->origin), SpeciesReference_getId(re[j]->products[l]->origin)) == 0){
								flag = 0;
								break;
							}
						}
						if(!flag){
							break;
						}
						for(l=0; l<re[j]->num_of_reactants; l++){
							if(SpeciesReference_isSetId(re[j]->reactants[l]->origin)
							   && strcmp(Rule_getVariable(rule[i]->origin), SpeciesReference_getId(re[j]->reactants[l]->origin)) == 0){
								flag = 0;
								break;
							}
						}
						if(!flag){
							break;
						}
					}
				}
				node = (ASTNode_t*)Rule_getMath(rule[i]->origin);
				node = ASTNode_deepCopy(node);
				TRACE(("original math : "));
				check_AST(node, NULL);
				alter_tree_structure(m, &node, NULL, 0, cp_AST);
				/* unit */
				if(rule[i]->target_species != NULL){
					if(rule[i]->target_species->is_amount
					   && !rule[i]->target_species->has_only_substance_units
					   && Compartment_getSpatialDimensions(rule[i]->target_species->locating_compartment->origin) != 0){
						assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(rule[i]->target_species->origin))), 0);
					}else if(rule[i]->target_species->is_concentration
							 && rule[i]->target_species->has_only_substance_units){
						assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(rule[i]->target_species->origin))), 1);
					}
				}
				/* unit */
				TRACE(("altered math : "));
				check_AST(node, NULL);
				math_length = connect_delayval_with_eq(m, rule[i]->eq, sp, param, comp, re, node, 0);
				/* ASTNode_free(node); */
				TRACE(("math\n"));
				check_math(rule[i]->eq);
			}
		}

		/* event */
		for(i=0; i<num_of_events; i++){
			node = (ASTNode_t*)Trigger_getMath(Event_getTrigger(ev[i]->origin));
			node = ASTNode_deepCopy(node);
			check_AST(node, NULL);
			alter_tree_structure(m, &node, NULL, 0, cp_AST);
			math_length = connect_delayval_with_eq(m, ev[i]->eq, sp, param, comp, re, node, 0);
			check_AST(node, NULL);
			if (Event_getDelay(ev[i]->origin) != NULL){
				node = (ASTNode_t*)Delay_getMath(ev[i]->event_delay->origin);
				node = ASTNode_deepCopy(node);
				check_AST(node, NULL);
				alter_tree_structure(m, &node, NULL, 0, cp_AST);
				math_length = connect_delayval_with_eq(m, ev[i]->event_delay->eq, sp, param, comp, re, node, 0);
			}
			for(j=0; j<Event_getNumEventAssignments(ev[i]->origin); j++){
				node = (ASTNode_t*)EventAssignment_getMath(ev[i]->assignments[j]->origin);
				node = ASTNode_deepCopy(node);
				check_AST(node, NULL);
				if(Event_getDelay(ev[i]->origin) != NULL && Event_getUseValuesFromTriggerTime(ev[i]->origin)){
					pre_ev_alter_tree_structure(&node, NULL, 0, (ASTNode_t*)Delay_getMath(Event_getDelay(ev[i]->origin)));
					TRACE(("after pre ev alter math : "));
					check_AST(node, NULL);
					ev_alter_tree_structure(m, &node, NULL, 0, cp_AST);
					TRACE(("after ev alter math : "));
					check_AST(node, NULL);
					post_ev_alter_tree_structure(m, &node, NULL, 0);
					TRACE(("after post ev alter math : "));
					check_AST(node, NULL);
				}else{
					alter_tree_structure(m, &node, NULL, 0, cp_AST);
				}
				/* unit */
				if(ev[i]->assignments[j]->target_species != NULL){
					if(ev[i]->assignments[j]->target_species->is_amount
					   && !ev[i]->assignments[j]->target_species->has_only_substance_units
					   && Compartment_getSpatialDimensions(ev[i]->assignments[j]->target_species->locating_compartment->origin) != 0){
						assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(ev[i]->assignments[j]->target_species->origin))), 0);
					}else if(ev[i]->assignments[j]->target_species->is_concentration
							 && ev[i]->assignments[j]->target_species->has_only_substance_units){
						assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(ev[i]->assignments[j]->target_species->origin))), 1);
					}
				}
				/* unit */
				check_AST(node, NULL);
				math_length = connect_delayval_with_eq(m, ev[i]->assignments[j]->eq, sp, param, comp, re, node, 0);
			}
		}
	}
	/* 2) reallocate event firing times
	   (note that new max index of array is smaller than new_max_index by 1)*/
	for(i=0; i<num_of_events; i++) {
		if(ev[i] != NULL && Event_getDelay(ev[i]->origin) != NULL){
			ev[i]->firing_times = (double*)realloc(ev[i]->firing_times, sizeof(double)*(int)new_max_index-1);
			for(k=max_index; k<new_max_index-1; k++){
				ev[i]->firing_times[k] = sim_time+1;
			}
		}
	}
}

void reallocate_result_objects(myResult* res, double** value_time_p_fordelay, unsigned  int max_result_index){
	double* head_time_e;

	/* save the value of max_index in Result */
	res->num_of_delay_rows = max_result_index * 2;

	/* reallocate value arrays in Result */
	head_time_e = (double *)realloc(res->values_time_fordelay, sizeof(double) * res->num_of_delay_rows);

	/* error check */
	if (!head_time_e) {
		fprintf(stderr, "failed in reallocation of time memory.\n");
		exit(1);
	}else{
		res->values_time_fordelay = head_time_e;
		head_time_e += max_result_index;
		*value_time_p_fordelay = head_time_e;
	}
}

unsigned int connect_delayval_with_eq(Model_t *m, equation *eq, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], ASTNode_t *node, int index){
	unsigned int i,j;
	int flag;
	const char* name;
	ASTNode_t *left, *right, *comp_node;
	if((ASTNode_getType(node) == AST_LOGICAL_AND
        || ASTNode_getType(node) == AST_LOGICAL_OR
        || ASTNode_getType(node) == AST_LOGICAL_XOR)
	   && ASTNode_getNumChildren(node) > 2){
		ASTNode_reduceToBinary(node);
	}
	if(ASTNode_getType(node) == AST_FUNCTION_DELAY){
		left = ASTNode_getLeftChild(node);
		comp_node = NULL;
		if(ASTNode_getType(left) != AST_NAME){
			comp_node = ASTNode_getRightChild(left);
			left = ASTNode_getLeftChild(left);
		}
		name = ASTNode_getName(left);
		flag = 1;
		/* species */
		for(i=0; i<Model_getNumSpecies(m); i++){
			if(strcmp(name, Species_getId(sp[i]->origin)) == 0){
				/* connection */
				eq->delay_number[index] = sp[i]->delay_val;
				if(comp_node != NULL){
					for(j=0; j<Model_getNumCompartments(m); j++){
						if(strcmp(ASTNode_getName(comp_node), Compartment_getId(comp[j]->origin)) == 0){
							/* connection */
							eq->delay_comp_size[index] = comp[j]->delay_val;
							break;
						}
					}
				}
				index++;
				flag = 0;
				break;
			}
		}
		if(flag){
			/* parameter */
			for(i=0; i<Model_getNumParameters(m); i++){
				if(strcmp(name, Parameter_getId(param[i]->origin)) == 0){
					/* connection */
					eq->delay_number[index] = param[i]->delay_val;
					index++;
					flag = 0;
					break;
				}
			}
		}
		if(flag){
			/* compartment */
			for(i=0; i<Model_getNumCompartments(m); i++){
				if(strcmp(name, Compartment_getId(comp[i]->origin)) == 0){
					/* connection */
					eq->delay_number[index] = comp[i]->delay_val;
					index++;
					flag = 0;
					break;
				}
			}
		}
		if(flag){
			for(i=0; i<Model_getNumReactions(m); i++){
				/* products */
				for(j=0; j<re[i]->num_of_products; j++){
					if(SpeciesReference_isSetId(re[i]->products[j]->origin)
					   && strcmp(name, SpeciesReference_getId(re[i]->products[j]->origin)) == 0){
						/* connection */
						eq->delay_number[index] = re[i]->products[j]->delay_val;
						index++;
						flag = 0;
						break;
					}
				}
				if(!flag){
					break;
				}
				/* reactants */
				for(j=0; j<re[i]->num_of_reactants; j++){
					if(SpeciesReference_isSetId(re[i]->reactants[j]->origin)
					   && strcmp(name, SpeciesReference_getId(re[i]->reactants[j]->origin)) == 0){
						/* connection */
						eq->delay_number[index] = re[i]->reactants[j]->delay_val;
						index++;
						flag = 0;
						break;
					}
				}
				if(!flag){
					break;
				}
			}
		}
	}else if((left=ASTNode_getLeftChild(node)) != NULL){
		index = connect_delayval_with_eq(m, eq, sp, param, comp, re, left, index);
	}
	if((right=ASTNode_getRightChild(node)) != NULL){
		index = connect_delayval_with_eq(m, eq, sp, param, comp, re, right, index);
	}

  if(ASTNode_isOperator(node)
      || ASTNode_isFunction(node)
      || ASTNode_isBoolean(node)){
    index++;
  }else if(ASTNode_getType(node) == AST_NAME){
    name = ASTNode_getName(node);
    flag = 1;
    for(i=0; i<Model_getNumSpecies(m); i++){
      if(strcmp(name, Species_getId(sp[i]->origin)) == 0){
        index++;
        flag = 0;
        break;
      }
    }
    if(flag){
      for(i=0; i<Model_getNumParameters(m); i++){
        if(strcmp(name, Parameter_getId(param[i]->origin)) == 0){
          index++;
          flag = 0;
          break;
        }
      }
    }
    if(flag){
      for(i=0; i<Model_getNumCompartments(m); i++){
        if(strcmp(name, Compartment_getId(comp[i]->origin)) == 0){
          index++;
          flag = 0;
          break;
        }
      }
    }
    if(flag){
      for(i=0; i<Model_getNumReactions(m); i++){
        for(j=0; j<re[i]->num_of_products; j++){
          if(SpeciesReference_isSetId(re[i]->products[j]->origin)
              && strcmp(name, SpeciesReference_getId(re[i]->products[j]->origin)) == 0){
            index++;
            flag = 0;
            break;
          }
        }
        if(!flag){
          break;
        }
        for(j=0; j<re[i]->num_of_reactants; j++){
          if(SpeciesReference_isSetId(re[i]->reactants[j]->origin)
              && strcmp(name, SpeciesReference_getId(re[i]->reactants[j]->origin)) == 0){
            index++;
            flag = 0;
            break;
          }
        }
        if(!flag){
          break;
        }
      }
    }
    if(flag){
      if(strcmp(name, "time") == 0
          || strcmp(name, "t") == 0
          || strcmp(name, "s") == 0){
        index++;
      }
    }
  }else if(ASTNode_getType(node) == AST_NAME_TIME){
    index++;
  }else if(ASTNode_getType(node) == AST_NAME_AVOGADRO){
    index++;
  }else if(ASTNode_getType(node) == AST_CONSTANT_E
      || ASTNode_getType(node) == AST_CONSTANT_PI){
    index++;
  }else{
    index++;
  }
	return index;
}
