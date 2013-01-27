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

myResult* bifurcation_analysis(Model_t *m, double sim_time, double dt, int print_interval, double time, int order, int print_amount, int use_lazy_method, int is_explicit, unsigned int num_of_species, unsigned int num_of_parameters, unsigned int num_of_compartments, unsigned int num_of_reactions, unsigned int num_of_rules, unsigned int num_of_events, unsigned int num_of_initialAssignments, mySpecies* mySp[], myParameter* myParam[], myCompartment* myComp[], myReaction* myRe[], myRule* myRu[], myEvent* myEv[], myInitialAssignment* myInitAssign[], myAlgebraicEquations* myAlgEq, timeVariantAssignments* timeVarAssign, allocated_memory* mem, copied_AST* cp_AST, myResult* result, myResult* rtn, boolean bif_param_is_local, char* sta_var_id, char* bif_param_id, double bif_param_min, double bif_param_max, double bif_param_stepsize, double transition_time) {
	unsigned int i, a, b;
	int j;
	/* max(min) of steady state */
	double local_max = 0.0;
	double local_min = 0.0;

	/* when bifurcation parameter is local*/
	double bif_param_value = 0.0;

	int sta_var_column = 0;
	int bif_param_column = 0;


	double init_max = 0.0;

	/* print BA result to file*/
	FILE* BAfp = NULL;

	/*for randam numbers*/
	dsfmt_t d;
	dsfmt_init_gen_rand(&d, (unsigned int) my_time(NULL));


	/* calculate the maximum of the state variable (to specify the range of initial value[range where the randam number is])*/
	for (j = 0; j < rtn->num_of_columns_sp; j++) {
		if (strcmp(sta_var_id, rtn->column_name_sp[j]) == 0) {
			sta_var_column = j;
			init_max = search_max(rtn, sta_var_column);
		}
	}
	/* substitute the bifurcation parameter value */
	for (i = 0; i < num_of_parameters; i++){
		if (strcmp(bif_param_id, Parameter_getId(myParam[i]->origin)) == 0){
			myParam[i]->value = bif_param_min;
			bif_param_column = i;
		}
	}
	if (bif_param_is_local) {
		for(a = 0; a < Model_getNumReactions(m); a++) {
			for(b = 0; b < ListOf_size((ListOf_t*)KineticLaw_getListOfParameters(Reaction_getKineticLaw(myRe[a]->origin))); b++){
				if (strcmp(bif_param_id, Parameter_getId(KineticLaw_getParameter(Reaction_getKineticLaw(myRe[a]->origin), b))) == 0) {
					free_mySBML_objects(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, mem, cp_AST);
					free_myResult(result);
					/* mySpecies *mySp[num_of_species]; */
					mySp = (mySpecies**)malloc(sizeof(mySpecies*) * num_of_species);
					/* myParameter *myParam[num_of_parameters]; */
					myParam = (myParameter**)malloc(sizeof(myParameter*) * num_of_parameters);
					/* myCompartment *myComp[num_of_compartments]; */
					myComp = (myCompartment**)malloc(sizeof(mySpecies*) * num_of_compartments);
					/* myReaction *myRe[num_of_reactions]; */
					myRe = (myReaction**)malloc(sizeof(myReaction*) * num_of_reactions);
					/* myRule *myRu[num_of_rules]; */
					myRu = (myRule**)malloc(sizeof(myRule*) * num_of_rules);
					/* myEvent *myEv[num_of_events]; */
					myEv = (myEvent**)malloc(sizeof(myEvent*) * num_of_events);
					/* myInitialAssignment *myInitAssign[num_of_initialAssignments]; */
					myInitAssign = (myInitialAssignment**)malloc(sizeof(myInitialAssignment*) * num_of_initialAssignments);
					mem = (allocated_memory*)malloc(sizeof(allocated_memory));
					mem->num_of_allocated_memory = 0;
					cp_AST = (copied_AST*)malloc(sizeof(copied_AST));
					cp_AST->num_of_copied_AST = 0;
					bif_param_value = bif_param_min;
					create_mySBML_objects_forBA(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, &myAlgEq, &timeVarAssign, sim_time, dt, &time, mem, cp_AST, bif_param_id, bif_param_value);
					break;
				}
			}
		}
	}
	if (!bif_param_is_local) {
		if ((BAfp = my_fopen(BAfp, "bifurcation_analysis.csv", "w")) != NULL) {
			while (myParam[bif_param_column]->value <= bif_param_max) {
				time = 0.0;
				mySp[sta_var_column]->value = init_max * dsfmt_genrand_close_open(&d);
				if (is_explicit == 1) {
					rtn = simulate_explicit(m, result, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval, &time, order, print_amount, mem);
				}else{
					rtn = simulate_implicit(m, result, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval, &time, order, use_lazy_method, print_amount, mem);
				}
				local_max = search_local_max(rtn, sta_var_column, transition_time, sim_time);
				local_min = search_local_min(rtn, sta_var_column, transition_time, sim_time);
				if(bif_param_value == bif_param_min) {
					fprintf(BAfp, "%s,%s,%s\n", bif_param_id, sta_var_id, sta_var_id);
				}
				fprintf(BAfp, "%.16g,%.16g,%.16g\n", bif_param_value, local_max, local_min);
				myParam[bif_param_column]->value += bif_param_stepsize;
			}
		}
	}else{
		if ((BAfp = my_fopen(BAfp, "bifurcation_analysis.csv", "w")) != NULL) {
			while (bif_param_value <= bif_param_max) {
				result = create_myResult(m, mySp, myParam, myComp, sim_time, dt, print_interval);
				time = 0.0;
				mySp[sta_var_column]->value = init_max * dsfmt_genrand_close_open(&d);
				if (is_explicit == 1) {
					rtn = simulate_explicit(m, result, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval, &time, order, print_amount, mem);
				  }else{
					rtn = simulate_implicit(m, result, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval, &time, order, use_lazy_method, print_amount, mem);
				}
				local_max = search_local_max(rtn, sta_var_column, transition_time, sim_time);
				local_min = search_local_min(rtn, sta_var_column, transition_time, sim_time);
				if(bif_param_value == bif_param_min) {
					fprintf(BAfp, "%s,%s,%s\n", bif_param_id, sta_var_id, sta_var_id);
				  }
				fprintf(BAfp, "%.16g,%.16g,%.16g\n", bif_param_value, local_max, local_min);
				free_mySBML_objects(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, mem, cp_AST);
				if (bif_param_value + bif_param_stepsize < bif_param_max) {
					free_myResult(result);
				}
				/* mySpecies *mySp[num_of_species]; */
				mySp = (mySpecies**)malloc(sizeof(mySpecies*) * num_of_species);
				/* myParameter *myParam[num_of_parameters]; */
				myParam = (myParameter**)malloc(sizeof(myParameter*) * num_of_parameters);
				/* myCompartment *myComp[num_of_compartments]; */
				myComp = (myCompartment**)malloc(sizeof(mySpecies*) * num_of_compartments);
				/* myReaction *myRe[num_of_reactions]; */
				myRe = (myReaction**)malloc(sizeof(myReaction*) * num_of_reactions);
				/* myRule *myRu[num_of_rules]; */
				myRu = (myRule**)malloc(sizeof(myRule*) * num_of_rules);
				/* myEvent *myEv[num_of_events]; */
				myEv = (myEvent**)malloc(sizeof(myEvent*) * num_of_events);
				/* myInitialAssignment *myInitAssign[num_of_initialAssignments]; */
				myInitAssign = (myInitialAssignment**)malloc(sizeof(myInitialAssignment*) * num_of_initialAssignments);
				mem = (allocated_memory*)malloc(sizeof(allocated_memory));
				mem->num_of_allocated_memory = 0;
				cp_AST = (copied_AST*)malloc(sizeof(copied_AST));
				cp_AST->num_of_copied_AST = 0;
				bif_param_value += bif_param_stepsize;
				create_mySBML_objects_forBA(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, &myAlgEq, &timeVarAssign, sim_time, dt, &time, mem, cp_AST, bif_param_id, bif_param_value);
				if (bif_param_value > bif_param_max) {
					free_mySBML_objects(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, myAlgEq, timeVarAssign, mem, cp_AST);
				}
			}
		}
	}
	fclose(BAfp);
	return rtn;
}
