#ifndef LibSBMLSim_h
#define LibSBMLSim_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <sbml/SBMLTypes.h>

#include "common.h"
#include "methods.h"
#include "version.h"

/* Boolean */
typedef enum _boolean { false, true } boolean;

/* structures for efficient simulation */
typedef struct _equation equation;
typedef struct _mySpecies mySpecies;
typedef struct _mySpeciesReference mySpeciesReference;
typedef struct _myParameter myParameter;
typedef struct _myCompartment myCompartment;
typedef struct _myReaction myReaction;
typedef struct _myRule myRule;
typedef struct _timeVariantAssignments timeVariantAssignments;
typedef struct _myEvent myEvent;
typedef struct _myEventAssignment myEventAssignment;
typedef struct _myDelay myDelay;
typedef struct _myInitialAssignment myInitialAssignment;
typedef struct _myASTNode myASTNode;
typedef struct _myAlgebraicEquations myAlgebraicEquations;
typedef struct _myAlgTargetSp myAlgTargetSp;
typedef struct _myAlgTargetParam myAlgTargetParam;
typedef struct _myAlgTargetComp myAlgTargetComp;
typedef struct _allocated_memory allocated_memory;
typedef struct _copied_AST copied_AST;
typedef struct _myResult myResult;

struct _equation{
  double *number[MAX_MATH_LENGTH];
  int operator[MAX_MATH_LENGTH];
  double **delay_number[MAX_MATH_LENGTH];
  double **delay_comp_size[MAX_MATH_LENGTH];
  equation *explicit_delay_eq[MAX_MATH_LENGTH];
  int math_length;
};

struct _mySpecies{
  Species_t *origin;
  double value;
  double temp_value;
  boolean is_amount;
  boolean is_concentration;
  int has_only_substance_units;
  myCompartment *locating_compartment;
  double k[4]; /* for runge kutta */
  double **delay_val;
  myRule *depending_rule;
  double prev_val[3]; /* previous values for multistep solution */
  double prev_k[3]; /* previous values for multistep solution */
};

struct _mySpeciesReference{
  mySpecies *mySp;
  SpeciesReference_t *origin;
  equation *eq; /* for l2v4 */
  double value;
  double temp_value;
  double k[4]; /* for runge kutta */
  double **delay_val;
  myRule *depending_rule;
  double prev_val[3];
  double prev_k[3];
};

struct _myParameter{
  Parameter_t* origin;
  double value;
  double temp_value;
  double k[4]; /* for runge kutta */
  double **delay_val;
  myRule *depending_rule;
  double prev_val[3]; /* previous values for multistep solution */
  double prev_k[3]; /* previous values for multistep solution */
};

struct _myCompartment{
  Compartment_t* origin;
  double value; /* compartment "size" value */
  double temp_value;
  double k[4]; /* for runge kutta */
  double **delay_val;
  myRule *depending_rule;
  double prev_val[3]; /* previous values for multistep solution */
  double prev_k[3]; /* previous values for multistep solution */
  mySpecies *including_species[MAX_INCLUDING_SPECIES];
  int num_of_including_species;
};

struct _myReaction{
  Reaction_t *origin;
  equation *eq;
  mySpeciesReference **products;
  int num_of_products;
  mySpeciesReference **reactants;
  int num_of_reactants;
  boolean is_fast;
  boolean is_reversible;
  equation *products_equili_numerator;
  equation *reactants_equili_numerator;
};

struct _myRule{
  Rule_t *origin;
  equation *eq;
  mySpecies *target_species;
  myParameter *target_parameter;
  myCompartment *target_compartment;
  mySpeciesReference *target_species_reference;
  boolean is_rate;
  boolean is_assignment;
  boolean is_algebraic;
};

struct _timeVariantAssignments{
  int num_of_time_variant_assignments;
  equation *eq[MAX_TIME_VARIANT_ASSIGNMENT];
  char *target_id[MAX_TIME_VARIANT_ASSIGNMENT];
};

struct _myEvent{
  Event_t *origin;
  equation *eq; /* condition equation */
  myEventAssignment** assignments;
  boolean is_able_to_fire;
  myDelay *event_delay;
  double *firing_times;
  int num_of_delayed_events_que;
  int next_firing_index;
  boolean is_persistent;
  equation *priority_eq;
};

struct _myEventAssignment{
  EventAssignment_t *origin;
  equation *eq; /* assignment equation */
  mySpecies *target_species;
  myParameter *target_parameter;
  myCompartment *target_compartment;
  mySpeciesReference *target_species_reference;
};

struct _myDelay{
  Delay_t *origin;
  equation *eq;
};

struct _myInitialAssignment{
  InitialAssignment_t *origin;
  equation *eq;
  mySpecies *target_species;
  myParameter *target_parameter;
  myCompartment *target_compartment;
  mySpeciesReference *target_species_reference;
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
  int num_of_algebraic_rules;
  int num_of_algebraic_variables;
  char *variables_id[MAX_ALGEBRAIC_VARIABLES];
  equation ***coefficient_matrix; /* use num_of_algebraic_equations > 1 */
  equation **constant_vector; /* use num_of_algebraic_equations > 1 */
  equation *coefficient; /* use num_of_algebraic_equations == 1 */
  equation *constant; /* use num_of_algebraic_equations == 1 */
  myAlgTargetSp *alg_target_species[MAX_ALGEBRAIC_VARIABLES];
  myAlgTargetParam *alg_target_parameter[MAX_ALGEBRAIC_VARIABLES];
  myAlgTargetComp *alg_target_compartment[MAX_ALGEBRAIC_VARIABLES];
  int num_of_alg_target_sp;
  int num_of_alg_target_param;
  int num_of_alg_target_comp;
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

struct _allocated_memory{
  double *memory[MAX_ALLOCATED_MEMORY];
  int num_of_allocated_memory;
};

struct _copied_AST{
  ASTNode_t *ast[MAX_COPIED_AST];
  int num_of_copied_AST;
};

struct _myResult{
  int num_of_rows;
  int num_of_columns_sp;
  int num_of_columns_param;
  int num_of_columns_comp;
  const char *column_name_time;
  const char **column_name_sp;
  const char **column_name_param;
  const char **column_name_comp;
  double *values_time;
  double *values_sp;
  double *values_param;
  double *values_comp;
};

/* Substitute kineticlaw local parameter node to simple Real value node in AST Tree
 * for efficient calculation) */
void set_local_para_as_value(ASTNode_t *node, KineticLaw_t *kl);

/* Get mathematical equations for calculation in reverse polish Notation */
int get_equation(Model_t *m, equation *eq, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], ASTNode_t *node, int index, double sim_time, double dt, double *time, myInitialAssignment *initAssign[], char *time_variant_target_id[], int num_of_time_variant_targets, timeVariantAssignments *timeVarAssign, allocated_memory *mem);

/* Calculate equations wrriten in reverse polish notation */
double calc(equation *eq, double dt, int cycle, double *reverse_time, int rk_order);

/* Calculate event equations wrriten in reverse polish notation */
void calc_event(myEvent *event[], int num_of_events, double dt, double time, int cycle, double *reverse_time);

/* numerical integration by explicit method(Runge Kutta and Adams-Bashforth) */
myResult* simulate_explicit(Model_t *m, myResult *result, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], myRule *rule[], myEvent *event[], myInitialAssignment *initAssign[], myAlgebraicEquations *algEq, timeVariantAssignments *timeVarAssign, double sim_time, double dt, int print_interval, double *time, int order, int print_amount, allocated_memory *mem);

/* LU decomposition */
int lu_decomposition(double **A, int *p, int N);

/* forward & backward substitution using LU matrix */
int lu_solve(double **A, int *p, int N, double *b);

/* numerical integration by implicit method(Adams-Moulton and Backward-Difference) */
myResult* simulate_implicit(Model_t *m, myResult *result, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], myRule *rule[], myEvent *event[], myInitialAssignment *initAssign[], myAlgebraicEquations *algEq, timeVariantAssignments *timeVarAssign, double sim_time, double dt, int print_interval, double *time, int order, int use_lazy_method, int print_amount, allocated_memory *mem);

/* get end_cycle */
int get_end_cycle(double sim_time, double dt);

/* create contents of myResult object */
void create_myResult_content(Model_t *m, myResult* result, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], double sim_time, double dt, int print_interval);

/* create my SBML obejects for efficient simulations */
void create_mySBML_objects(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], myReaction *myRe[], myRule *myRu[], myEvent *myEv[], myInitialAssignment *myInitAssign[], myAlgebraicEquations **myAlgEq, timeVariantAssignments **timeVarAssign, double sim_time, double dt, double *time, allocated_memory *mem, copied_AST *cp_AST);

/* free all created my SBML objects */
void free_mySBML_objects(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], myReaction *myRe[], myRule *myRu[], myEvent *myEv[], myInitialAssignment *myInitAssign[], myAlgebraicEquations *myAlgEq, timeVariantAssignments *timeVarAssign, double sim_time, double dt, allocated_memory *mem, copied_AST *cp_AST);

/* print result column list */
void print_result_list(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[]);

/* print result */
void print_result(myResult* result);
void write_result(myResult* result, char* file);
void write_csv(myResult* result, char* file);
void print_result_to_file(myResult* result, char* file, char delimiter);
void output_result(myResult* result, FILE* fp, char delimiter);
void write_separate_result(myResult* result, char* file_s, char* file_p, char* file_c); 

/* calc k(gradient or value of algebraic or assignment rule) */
void calc_k(mySpecies *sp[], int sp_num, myParameter *param[], int param_num, myCompartment *comp[], int comp_num, mySpeciesReference *spr[], int spr_num, myReaction *re[], int re_num, myRule *rule[], int rule_num, int cycle, double dt, double *reverse_time, int use_rk, int call_first_time_in_cycle);

/* calculate temp_value using k */
void calc_temp_value(mySpecies *sp[], int sp_num, myParameter *param[], int param_num, myCompartment *comp[], int comp_num, mySpeciesReference *spr[], int spr_num, int cycle, double dt, int use_rk);

/* forwarding value(substitute calculated temp value to value) */
void forwarding_value(mySpecies *sp[], int sp_num, myParameter *param[], int param_num, myCompartment *comp[], int comp_num, mySpeciesReference *spr[], int spr_num);

/* check the number of object for culculation */
void check_num(int num_of_species, int num_of_parameters, int num_of_compartments, int num_of_reactions, int *num_of_all_var_species, int *num_of_all_var_parameters, int *num_of_all_var_compartments, int *num_of_all_var_species_reference, int *num_of_var_species, int *num_of_var_parameters, int *num_of_var_compartments, int *num_of_var_species_reference, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[]);

/* create list of object for calculation */
void create_calc_object_list(int num_of_species, int num_of_parameters, int num_of_compartments, int num_of_reactions, int num_of_all_var_species, int num_of_all_var_parameters, int num_of_all_var_compartments, int num_of_all_var_species_reference, int num_of_var_species, int num_of_var_parameters, int num_of_var_compartments, int num_of_var_species_reference, mySpecies *all_var_sp[], myParameter *all_var_param[], myCompartment *all_var_comp[], mySpeciesReference *all_var_spr[], mySpecies *var_sp[], myParameter *var_param[], myCompartment *var_comp[], mySpeciesReference *var_spr[], mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[]);

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
void prepare_algebraic(Model_t *m, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], myRule *ru[], myEvent *ev[], myInitialAssignment *initAssign[], myAlgebraicEquations *algEq, double sim_time, double dt, double *time, char *time_variant_target_id[], int num_of_time_variant_targets, timeVariantAssignments *timeVarAssign, allocated_memory *mem, copied_AST *cp_AST);

void alg_alter_tree_structure(ASTNode_t **node_p, ASTNode_t *parent, int child_order);

/* event alter */
void pre_ev_alter_tree_structure(ASTNode_t **node_p, ASTNode_t *parent, int child_order, ASTNode_t *delay_math);
void ev_alter_tree_structure(Model_t *m, ASTNode_t **node_p, ASTNode_t *parent, int child_order, copied_AST *cp_AST);
void post_ev_alter_tree_structure(Model_t *m, ASTNode_t **node_p, ASTNode_t *parent, int child_order);

/* assignment alter */
void assignment_alter_tree_structure(ASTNode_t **node_p, char* comp_name, int sw);

/* myASTNode_func */
void myASTNode_create(myASTNode *myNode, ASTNode_t *node, myASTNode *copied_myAST[], int *num_of_copied_myAST);

void ASTNode_recreate(myASTNode *myNode, ASTNode_t *node);

void myASTNode_free(myASTNode *copied_myAST[], int num_of_copied_myAST);

void check_myAST(myASTNode *myNode);

/* prepare reversible fast reaction */
void prepare_reversible_fast_reaction(Model_t *m, myReaction *re[], mySpecies *sp[], myParameter *param[], myCompartment *comp[], double sim_time, double dt, double *time, myInitialAssignment *initAssign[], char *time_variant_target_id[], int num_of_time_variant_targets, timeVariantAssignments *timeVarAssign, allocated_memory *mem, copied_AST *cp_AST);

/* calculate initial assignment */
void calc_initial_assignment(myInitialAssignment *initAssign[], int num_of_initialAssignments, double dt, int cycle, double *reverse_time);

/* initialize delay val */
void initialize_delay_val(mySpecies *sp[], int num_of_species, myParameter *param[], int num_of_parameters, myCompartment *comp[], int num_of_compartments, myReaction *re[], int num_of_reactions, double sim_time, double dt, int last_call);

/* substitute delay val */
void substitute_delay_val(mySpecies *sp[], int num_of_species, myParameter *param[], int num_of_parameters, myCompartment *comp[], int num_of_compartments, myReaction *re[], int num_of_reactions, int cycle);

/* Run Simulation and output CSV */
myResult* simulateSBMLModel(Model_t *m, myResult *result, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method);

/* Run Simulation from SBML string */
myResult* simulateSBMLFromString(const char* str, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method);

/* following functions are defined in math_functions.c */
int64_t factorial(int n);

#endif  /* LibSBMLSim_h */
