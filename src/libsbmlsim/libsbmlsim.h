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
#ifndef LibSBMLSim_h
#define LibSBMLSim_h

#ifdef _MSC_VER
#define _CRT_SECURE_NO_DEPRECATE
#endif /* _MSC_VER */

#include <sbml/SBMLTypes.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif /* _MSC_VER */
#include <math.h>
#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288
#endif /* M_PI */
#ifndef M_E
#define M_E 2.71828182845904523536028747135266250
#endif /* M_E */

#include <float.h>
#include <time.h>

#include "osarch.h"
#include "common.h"
#include "methods.h"
#include "errorcodes.h"
#include "boolean.h"
#include "version.h"

#include "typedefs.h"
#include "equation.h"
#include "myResult.h"
#include "mySpecies.h"
#include "mySpeciesReference.h"
#include "myParameter.h"
#include "myCompartment.h"
#include "myReaction.h"
#include "myEvent.h"
#include "myEventAssignment.h"
#include "myInitialAssignment.h"
#include "myRule.h"
#include "myDelay.h"
#include "allocated_memory.h"
#include "copied_AST.h"

#define DSFMT_MEXP 19937
#include "dSFMT-params19937.h"
#include "dSFMT.h"
#include "dSFMT-params.h"


struct _timeVariantAssignments{
  unsigned int num_of_time_variant_assignments;
  equation *eq[MAX_TIME_VARIANT_ASSIGNMENT];
  char *target_id[MAX_TIME_VARIANT_ASSIGNMENT];
};

struct _myASTNode{
  ASTNode_t *origin;
  myASTNode *parent;
  myASTNode *left;
  myASTNode *right;
};

/* num_of_algebraic_equation and num_of_algebraic_variables
 * must be the same
 * coefficinet matrix size must be regular matrix and the size of this is
 * num_of_algebraic_variables * num_of_algebraic_variables)  */
struct _myAlgebraicEquations{
  unsigned int num_of_algebraic_rules;
  unsigned int num_of_algebraic_variables;
  char *variables_id[MAX_ALGEBRAIC_VARIABLES];
  equation ***coefficient_matrix; /* use num_of_algebraic_equations > 1 */
  equation **constant_vector; /* use num_of_algebraic_equations > 1 */
  equation *coefficient; /* use num_of_algebraic_equations == 1 */
  equation *constant; /* use num_of_algebraic_equations == 1 */
  myAlgTargetSp *alg_target_species[MAX_ALGEBRAIC_VARIABLES];
  myAlgTargetParam *alg_target_parameter[MAX_ALGEBRAIC_VARIABLES];
  myAlgTargetComp *alg_target_compartment[MAX_ALGEBRAIC_VARIABLES];
  unsigned int num_of_alg_target_sp;
  unsigned int num_of_alg_target_param;
  unsigned int num_of_alg_target_comp;
  mySpecies *target_species;
  myParameter *target_parameter;
  myCompartment *target_compartment;
};

struct _myAlgTargetSp{
  mySpecies *target_species;
  int order;
};

struct _myAlgTargetParam{
  myParameter *target_parameter;
  int order;
};

struct _myAlgTargetComp{
  myCompartment *target_compartment;
  int order;
};

/* Substitute kineticlaw local parameter node to simple Real value node in AST Tree
 * for efficient calculation) */
void set_local_para_as_value(ASTNode_t *node, KineticLaw_t *kl);
void set_local_para_as_value_forBA(ASTNode_t *node, KineticLaw_t *kl, char* bif_param_id, double bif_param_value);


/* Get mathematical equations for calculation in reverse polish Notation */
unsigned int get_equation(Model_t *m, equation *eq, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], ASTNode_t *node, unsigned int index, double sim_time, double dt, double *time, myInitialAssignment *initAssign[], char *time_variant_target_id[], unsigned int num_of_time_variant_targets, timeVariantAssignments *timeVarAssign, allocated_memory *mem);

unsigned int get_equationf(Model_t *m, equation *eq, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], ASTNode_t *node, unsigned int index, double sim_time, double dt, double *time, myInitialAssignment *initAssign[], char *time_variant_target_id[], unsigned int num_of_time_variant_targets, timeVariantAssignments *timeVarAssign, allocated_memory *mem, int print_interval);


/* Calculate equations wrriten in reverse polish notation */
double calc(equation *eq, double dt, int cycle, double *reverse_time, int rk_order);

double calcf(equation *eq, double dt, int cycle, double *reverse_time, int rk_order, double* time, double* stage_time, myResult* res, int print_interval, int* err_zero_flag);




/* Calculate event equations wrriten in reverse polish notation */
void calc_event(myEvent *event[], unsigned int num_of_events, double dt, double time, int cycle, double *reverse_time);

void calc_eventf(myEvent *event[], unsigned int num_of_events, double dt, double time, int cycle, double *reverse_time, myResult* res, int print_interval, int* err_zero_flag);

void recursive_calc_event(myEvent *event[], unsigned int num_of_events, myEvent *event_buf[], unsigned int *num_of_remained_events, double *assignment_values_from_trigger_time[], double dt, double time, int cycle, double *reverse_time);

void recursive_calc_eventf(myEvent *event[], unsigned int num_of_events, myEvent *event_buf[], unsigned int *num_of_remained_events, double *assignment_values_from_trigger_time[], double dt, double time, int cycle, double *reverse_time, myResult* res, int print_interval, int* err_zero_flag);

/* numerical integration by explicit method(Runge Kutta and Adams-Bashforth) */
myResult* simulate_explicit(Model_t *m, myResult *result, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], myRule *rule[], myEvent *event[], myInitialAssignment *initAssign[], myAlgebraicEquations *algEq, timeVariantAssignments *timeVarAssign, double sim_time, double dt, int print_interval, double *time, int order, int print_amount, allocated_memory *mem);

myResult* simulate_explicitf(Model_t *m, myResult* result, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], myRule *rule[], myEvent *event[], myInitialAssignment *initAssign[], myAlgebraicEquations *algEq, timeVariantAssignments *timeVarAssign, double sim_time, double dt, int print_interval, double *time, int order, int print_amount, allocated_memory *mem, double atol, double rtol, double facmax, copied_AST *cp_AST, int* err_zero_flag);

/* count the number of ODE [for variable stepsize] */
int count_ode(mySpecies* sp[], unsigned int num_of_species, int* ode_check, Species_t* s);



/* LU decomposition */
int lu_decomposition(double **A, int *p, int N);

/* forward & backward substitution using LU matrix */
int lu_solve(double **A, int *p, int N, double *b);

/* numerical integration by implicit method(Adams-Moulton and Backward-Difference) */
myResult* simulate_implicit(Model_t *m, myResult *result, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], myRule *rule[], myEvent *event[], myInitialAssignment *initAssign[], myAlgebraicEquations *algEq, timeVariantAssignments *timeVarAssign, double sim_time, double dt, int print_interval, double *time, int order, int use_lazy_method, int print_amount, allocated_memory *mem);

/** util.c **/
/* get end_cycle */
int get_end_cycle(double sim_time, double dt);

/* set seed for random */
void set_seed(void);

/* my_time function. We created this function because we can't use time() and variable time in the same file. Argh... */
time_t my_time(time_t* tloc);

/* strdup is not supported in C89, so I reimplement it. */
char* dupstr(const char *str);

/* create myResult object (and contents) */
myResult *create_myResult(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], double sim_time, double dt, int print_interval);

myResult *create_myResultf(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], double sim_time, double dt);

/* create myResult object with error. */
/* create_myResult_with_errorCode insert default error_message */
myResult *create_myResult_with_error(LibsbmlsimErrorCode code, const char *message);
myResult *create_myResult_with_errorCode(LibsbmlsimErrorCode code);

/* return 1 if error */
SBMLSIM_EXPORT int myResult_isError(myResult *result);

/* get error message */
SBMLSIM_EXPORT const char *myResult_getErrorMessage(myResult *result);

/* deallocate myResult */
SBMLSIM_EXPORT void free_myResult(myResult *res);
SBMLSIM_EXPORT void __free_myResult(myResult *res);

/* create my SBML obejects for efficient simulations */
void create_mySBML_objects(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], myReaction *myRe[], myRule *myRu[], myEvent *myEv[], myInitialAssignment *myInitAssign[], myAlgebraicEquations **myAlgEq, timeVariantAssignments **timeVarAssign, double sim_time, double dt, double *time, allocated_memory *mem, copied_AST *cp_AST);

void create_mySBML_objects_forBA(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], myReaction *myRe[], myRule *myRu[], myEvent *myEv[], myInitialAssignment *myInitAssign[], myAlgebraicEquations **myAlgEq, timeVariantAssignments **timeVarAssign, double sim_time, double dt, double *time, allocated_memory *mem, copied_AST *cp_AST, char* bif_param_id, double bif_param_value);

void create_mySBML_objectsf(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], myReaction *myRe[], myRule *myRu[], myEvent *myEv[], myInitialAssignment *myInitAssign[], myAlgebraicEquations **myAlgEq, timeVariantAssignments **timeVarAssign, double sim_time, double dt, double *time, allocated_memory *mem, copied_AST *cp_AST, int print_interval);


/* free all created my SBML objects */
void free_mySBML_objects(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], myReaction *myRe[], myRule *myRu[], myEvent *myEv[], myInitialAssignment *myInitAssign[], myAlgebraicEquations *myAlgEq, timeVariantAssignments *timeVarAssign, allocated_memory *mem, copied_AST *cp_AST);

/* print result column list */
void print_result_list(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[]);

/* show species (parameter) list [for bifurcation analysis]*/
void show_sp(Model_t *m);
void show_para(Model_t *m);

/* return the maximum of state variable [for bifurcation analysis]*/
 double search_max(myResult* result, int sta_var_column);
/* return the local maximum or minimum of state variable [for bifurcation analysis]*/
double search_local_max(myResult* result, int sta_var_column, double transition_time, double sim_time);
double search_local_min(myResult* result, int sta_var_column, double transition_time, double sim_time);
/* print result */
FILE* my_fopen(FILE* fp, const char* file, char* mode);
SBMLSIM_EXPORT void print_result(myResult* result);
SBMLSIM_EXPORT void write_result(myResult* result, const char* file);
SBMLSIM_EXPORT void write_csv(myResult* result, const char* file);
void print_result_to_file(myResult* result, const char* file, char delimiter);
void output_result(myResult* result, FILE* fp, char delimiter);
SBMLSIM_EXPORT void write_separate_result(myResult* result, const char* file_s, const char* file_p, const char* file_c); 

/* calc k(gradient or value of algebraic or assignment rule) */
void calc_k(mySpecies *sp[], unsigned int sp_num, myParameter *param[], unsigned int param_num, myCompartment *comp[], unsigned int comp_num, mySpeciesReference *spr[], unsigned int spr_num, myReaction *re[], unsigned int re_num, myRule *rule[], unsigned int rule_num, int cycle, double dt, double *reverse_time, int use_rk, int call_first_time_in_cycle);

void calc_kf(mySpecies *sp[], unsigned int sp_num, myParameter *param[], unsigned int param_num, myCompartment *comp[], unsigned int comp_num, mySpeciesReference *spr[], unsigned int spr_num, myReaction *re[], unsigned int re_num, myRule *rule[], unsigned int rule_num, int cycle, double dt, double *reverse_time, int use_rk, int call_first_time_in_cycle, double* time, myResult* res, myAlgebraicEquations *algEq, int print_interval, int* err_zero_flag, int order);

int calc_by_algebraic(myAlgebraicEquations *algEq, int cycle, double dt, double reverse_time, double* time, myResult* result, int print_interval, int* err_zero_flag);

void calc_by_assignment(mySpecies *sp[], unsigned int sp_num, myParameter *param[], unsigned int param_num, myCompartment *comp[], unsigned int comp_num, mySpeciesReference *spr[], unsigned int spr_num, double dt, int cycle, double reverse_time, double* time, myResult* result, int print_interval, int* err_zero_flag);


/* calculate temp_value using k */
void calc_temp_value(mySpecies *sp[], unsigned int sp_num, myParameter *param[], unsigned int param_num, myCompartment *comp[], unsigned int comp_num, mySpeciesReference *spr[], unsigned int spr_num, double dt, int use_rk);

double calc_sum_error(mySpecies *sp[], unsigned int sp_num, myParameter *param[], unsigned int param_num, myCompartment *comp[], unsigned int comp_num, mySpeciesReference *spr[], unsigned int spr_num, double dt, int use_rk, double atol, double rtol, int* ode_num, int time_progressed_flag, int order);

double calc_error(double dxdt, double dxdt4, double cur_value, double next_value, double atol, double rtol);

void free_dxdt(double *sp_dxdt, double *sp_dxdt4, double *param_dxdt, double *param_dxdt4, double *comp_dxdt, double *comp_dxdt4, double *spr_dxdt, double *spr_dxdt4);

double calc_eps(double value);

/* forwarding value(substitute calculated temp value to value) */
void forwarding_value(mySpecies *sp[], int sp_num, myParameter *param[], int param_num, myCompartment *comp[], int comp_num, mySpeciesReference *spr[], int spr_num);

/* check the number of object for culculation */
void check_num(unsigned int num_of_species, unsigned int num_of_parameters, unsigned int num_of_compartments, unsigned int num_of_reactions, unsigned int *num_of_all_var_species, unsigned int *num_of_all_var_parameters, unsigned int *num_of_all_var_compartments, unsigned int *num_of_all_var_species_reference, unsigned int *num_of_var_species, unsigned int *num_of_var_parameters, unsigned int *num_of_var_compartments, unsigned int *num_of_var_species_reference, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[]);

/* create list of object for calculation */
void create_calc_object_list(unsigned int num_of_species, unsigned int num_of_parameters, unsigned int num_of_compartments, unsigned int num_of_reactions, mySpecies *all_var_sp[], myParameter *all_var_param[], myCompartment *all_var_comp[], mySpeciesReference *all_var_spr[], mySpecies *var_sp[], myParameter *var_param[], myCompartment *var_comp[], mySpeciesReference *var_spr[], mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[]);

/* function for checking whether string is number */
boolean str_is_number(const char *str);

/* chomp */
void chomp(char *str);

/* Alter the AST structure for calculation */
void alter_tree_structure(Model_t *m, ASTNode_t **node_p, ASTNode_t *parent, int child_order, copied_AST *cp_AST);

/* Checker for reverse polish notation */
void check_math(equation *eq);

/* Checker for AST */
void check_AST(ASTNode_t *node, ASTNode_t *parent);

/* print node type */
void print_node_type(ASTNode_t *node);

/* function for algebraic preparation */
void prepare_algebraic(Model_t *m, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], myRule *ru[], myEvent *ev[], myInitialAssignment *initAssign[], myAlgebraicEquations *algEq, double sim_time, double dt, double *time, char *time_variant_target_id[], unsigned int num_of_time_variant_targets, timeVariantAssignments *timeVarAssign, allocated_memory *mem, copied_AST *cp_AST);

void alg_alter_tree_structure(ASTNode_t **node_p, ASTNode_t *parent, int child_order);

/* event alter */
void pre_ev_alter_tree_structure(ASTNode_t **node_p, ASTNode_t *parent, int child_order, ASTNode_t *delay_math);
void ev_alter_tree_structure(Model_t *m, ASTNode_t **node_p, ASTNode_t *parent, int child_order, copied_AST *cp_AST);
void post_ev_alter_tree_structure(Model_t *m, ASTNode_t **node_p, ASTNode_t *parent, int child_order);

/* assignment alter */
void assignment_alter_tree_structure(ASTNode_t **node_p, char* comp_name, int sw);

/* myASTNode_func */
void myASTNode_create(myASTNode *myNode, ASTNode_t *node, myASTNode *copied_myAST[], unsigned int *num_of_copied_myAST);

void ASTNode_recreate(myASTNode *myNode, ASTNode_t *node);

void myASTNode_free(myASTNode *copied_myAST[], unsigned int num_of_copied_myAST);

void check_myAST(myASTNode *myNode);

/* prepare reversible fast reaction */
void prepare_reversible_fast_reaction(Model_t *m, myReaction *re[], mySpecies *sp[], myParameter *param[], myCompartment *comp[], double sim_time, double dt, double *time, myInitialAssignment *initAssign[], char *time_variant_target_id[], unsigned int num_of_time_variant_targets, timeVariantAssignments *timeVarAssign, allocated_memory *mem, copied_AST *cp_AST);

/* calculate initial assignment */
void calc_initial_assignment(myInitialAssignment *initAssign[], unsigned int num_of_initialAssignments, double dt, int cycle, double *reverse_time);

void calc_initial_assignmentf(myInitialAssignment *initAssign[], unsigned int num_of_initialAssignments, double dt, int cycle, double *reverse_time, double* time, myResult* result, int print_interval, int* err_zero_flag);

/* initialize delay val */
void initialize_delay_val(mySpecies *sp[], unsigned int num_of_species, myParameter *param[], unsigned int num_of_parameters, myCompartment *comp[], unsigned int num_of_compartments, myReaction *re[], unsigned int num_of_reactions, double sim_time, double dt, int last_call);

void initialize_delay_valf(mySpecies *sp[], unsigned int num_of_species, myParameter *param[], unsigned int num_of_parameters, myCompartment *comp[], unsigned int num_of_compartments, myReaction *re[], unsigned int num_of_reactions, double sim_time, double dt, int last_call);


/* substitute delay val */
void substitute_delay_val(mySpecies *sp[], unsigned int num_of_species, myParameter *param[], unsigned int num_of_parameters, myCompartment *comp[], unsigned int num_of_compartments, myReaction *re[], unsigned int num_of_reactions, int cycle);

void substitute_delay_valf(mySpecies *sp[], unsigned int num_of_species, myParameter *param[], unsigned int num_of_parameters, myCompartment *comp[], unsigned int num_of_compartments, myReaction *re[], unsigned int num_of_reactions, int cycle);

/* Debug print */
void dbg_printf(const char *fmt, ...);

/* Progress print */
void prg_printf(const char *fmt, ...);

/* Run Simulation from SBML Model */
SBMLSIM_EXPORT myResult* simulateSBMLModel(Model_t *m, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method, double atol, double rtol, double facmax);

SBMLSIM_EXPORT myResult* simulateSBMLModelf(Model_t *m, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method, double atol, double rtol, double facmax, int use_bifurcation_analysis);

/* Run Simulation from SBML string */
SBMLSIM_EXPORT myResult* simulateSBMLFromString(const char* str, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method);

/* Run Simulation from SBML file */
SBMLSIM_EXPORT myResult* simulateSBMLFromFile(const char* file, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method);

/* Bifurcation Analysis mode */
myResult* bifurcation_analysis(Model_t *m, double sim_time, double dt, int print_interval, double time, int order, int print_amount, int use_lazy_method, int is_explicit, unsigned int num_of_species, unsigned int num_of_parameters, unsigned int num_of_compartments, unsigned int num_of_reactions, unsigned int num_of_rules, unsigned int num_of_events, unsigned int num_of_initialAssignments, mySpecies* mySp[], myParameter* myParam[], myCompartment* myComp[], myReaction* myRe[], myRule* myRu[], myEvent* myEv[], myInitialAssignment* myInitAssign[], myAlgebraicEquations* myAlgEq, timeVariantAssignments* timeVarAssign, allocated_memory* mem, copied_AST* cp_AST, myResult* result, myResult* rtn, boolean bif_param_is_local, char* sta_var_id, char* bif_param_id, double bif_param_min, double bif_param_max, double bif_param_stepsize, double transition_time);

/*for vaariable step-size integration */
/* calculate the solution in the past by linear approximation */
double approximate_delay_linearly(double* stack, int pos, double** delay_preserver, double* time, int rk_order, myResult* res, int cycle, int print_interval, int* err_zero_flag);

/* rearrange calculation result by linear approximation*/
double approximate_printresult_linearly(double value, double temp_value, double value_time, double tempvalue_time, double fixed_time);

/* reallocate objects (for variable step-size)*/
void realloc_mySBML_objects(Model_t *m, mySpecies *sp[], unsigned int num_of_species, myParameter *param[], unsigned int num_of_parameters, myCompartment *comp[], unsigned int num_of_compartments, myReaction *re[], unsigned int num_of_reactions, myRule *rule[], myEvent *ev[], unsigned int num_of_events, myInitialAssignment *myInitAssign[], timeVariantAssignments **timeVarAssign, copied_AST *cp_AST, double sim_time, int max_index);

unsigned int connect_delayval_with_eq(Model_t *m, equation *eq, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], ASTNode_t *node, int index);

void reallocate_result_objects(myResult* res, double** value_time_p_fordelay, unsigned  int max_result_index);


/** math_functions.c **/
int64_t factorial(int n);
double my_fmax(double a, double b);
double my_fmin(double a, double b);
double my_asinh(double x);
double my_acosh(double x);
double my_atanh(double x);
int my_isnan(double x);
#if defined(_MSC_VER) || defined(__STRICT_ANSI__)
double s_log1p(double x);
double __ieee754_acosh(double x);
double __ieee754_atanh(double x);
#endif

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  /* LibSBMLSim_h */
