#include "sbml/common/common.h"
#include "sbml/SBMLReader.h"
#include "sbml/SBMLTypes.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include "header.h"

int get_equation(Model_t *m, equation *eq, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], ASTNode_t *node, int index, double sim_time, double dt, double *time, myInitialAssignment *initAssign[], char *time_variant_target_id[], int num_of_time_variant_targets, timeVariantAssignments *timeVarAssign, allocated_memory *mem){
  int i, j, k, flag;
  int op;
  const char *name;
  double value;
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
    for(i=0; i<Model_getNumSpecies(m); i++){
      if(strcmp(name, Species_getId(sp[i]->origin)) == 0){
	//create delay
	if(sp[i]->delay_val == NULL){
	  sp[i]->delay_val = (double**)malloc(sizeof(double*)*(int)(sim_time/dt+1));
	  for(j=0; j<(int)(sim_time/dt+1); j++){
	    sp[i]->delay_val[j] = (double*)malloc(sizeof(double)*4);
	  }
	}
	if(comp_node != NULL){
	  dbg_printf("comp delay creation for species start\n");
	  for(j=0; j<Model_getNumCompartments(m); j++){
	    if(strcmp(ASTNode_getName(comp_node), Compartment_getId(comp[j]->origin)) == 0){
	      if(comp[j]->delay_val == NULL){
		comp[j]->delay_val = (double**)malloc(sizeof(double*)*(int)(sim_time/dt+1));
		for(k=0; k<(int)(sim_time/dt+1); k++){
		  comp[j]->delay_val[k] = (double*)malloc(sizeof(double)*4);
		}
	      }
	    }
	  }
	}
	dbg_printf("comp delay creation for species finish\n");
	eq->number[index] = NULL;
	eq->operator[index] = 0;
	eq->delay_number[index] = sp[i]->delay_val;
	if(comp_node != NULL){
	  for(j=0; j<Model_getNumCompartments(m); j++){
	    if(strcmp(ASTNode_getName(comp_node), Compartment_getId(comp[j]->origin)) == 0){
	      eq->delay_comp_size[index] = comp[j]->delay_val;
	      break;
	    }
	  }
	}else{
	  eq->delay_comp_size[index] = NULL;
	}
	eq->explicit_delay_eq[index] = NULL;
	if(initAssign != NULL){
	  for(j=0; j<num_of_time_variant_targets; j++){
	    if(strcmp(time_variant_target_id[j], name) == 0){
	      for(k=0; k<Model_getNumInitialAssignments(m); k++){
		if(strcmp(InitialAssignment_getSymbol(initAssign[k]->origin), name) == 0){
		  eq->explicit_delay_eq[index] = initAssign[k]->eq;
		}
	      }
	    }
	  }
	}
	if(timeVarAssign != NULL){
	  for(j=0; j<timeVarAssign->num_of_time_variant_assignments; j++){
	    if(strcmp(timeVarAssign->target_id[j], name) == 0){
	      eq->explicit_delay_eq[index] = timeVarAssign->eq[j];
	    }
	  }
	}
	index++;
	flag = 0;
	break;
      }
    }
    if(flag){
      for(i=0; i<Model_getNumParameters(m); i++){
	if(strcmp(name, Parameter_getId(param[i]->origin)) == 0){
	  //create delay
	  if(param[i]->delay_val == NULL){
	    param[i]->delay_val = (double**)malloc(sizeof(double*)*(int)(sim_time/dt+1));
	    for(j=0; j<(int)(sim_time/dt+1); j++){
	      param[i]->delay_val[j] = (double*)malloc(sizeof(double)*4);
	    }
	  }
	  eq->number[index] = NULL;
	  eq->operator[index] = 0;
	  eq->delay_number[index] = param[i]->delay_val;
	  eq->delay_comp_size[index] = NULL;
	  eq->explicit_delay_eq[index] = NULL;
	  if(initAssign != NULL){
	    for(j=0; j<num_of_time_variant_targets; j++){
	      if(strcmp(time_variant_target_id[j], name) == 0){
		for(k=0; k<Model_getNumInitialAssignments(m); k++){
		  if(strcmp(InitialAssignment_getSymbol(initAssign[k]->origin), name) == 0){
		    eq->explicit_delay_eq[index] = initAssign[k]->eq;
		  }
		}
	      }
	    }
	  }
	  if(timeVarAssign != NULL){
	    for(j=0; j<timeVarAssign->num_of_time_variant_assignments; j++){
	      if(strcmp(timeVarAssign->target_id[j], name) == 0){
		eq->explicit_delay_eq[index] = timeVarAssign->eq[j];
	      }
	    }
	  }
	  index++;
	  flag = 0;
	  break;
	}
      }
    }
    if(flag){
      for(i=0; i<Model_getNumCompartments(m); i++){
	if(strcmp(name, Compartment_getId(comp[i]->origin)) == 0){
	  //create delay
	  if(comp[i]->delay_val == NULL){
	    comp[i]->delay_val = (double**)malloc(sizeof(double*)*(int)(sim_time/dt+1));
	    for(j=0; j<(int)(sim_time/dt+1); j++){
	      comp[i]->delay_val[j] = (double*)malloc(sizeof(double)*4);
	    }
	  }
	  eq->number[index] = NULL;
	  eq->operator[index] = 0;
	  eq->delay_number[index] = comp[i]->delay_val;
	  eq->delay_comp_size[index] = NULL;
	  eq->explicit_delay_eq[index] = NULL;
	  if(initAssign != NULL){
	    for(j=0; j<num_of_time_variant_targets; j++){
	      if(strcmp(time_variant_target_id[j], name) == 0){
		for(k=0; k<Model_getNumInitialAssignments(m); k++){
		  if(strcmp(InitialAssignment_getSymbol(initAssign[k]->origin), name) == 0){
		    eq->explicit_delay_eq[index] = initAssign[k]->eq;
		  }
		}
	      }
	    }
	  }
	  if(timeVarAssign != NULL){
	    for(j=0; j<timeVarAssign->num_of_time_variant_assignments; j++){
	      if(strcmp(timeVarAssign->target_id[j], name) == 0){
		eq->explicit_delay_eq[index] = timeVarAssign->eq[j];
	      }
	    }
	  }
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
	    //create delay
	    if(re[i]->products[j]->delay_val == NULL){
	      re[i]->products[j]->delay_val = (double**)malloc(sizeof(double*)*(int)(sim_time/dt+1));
	      for(k=0; k<(int)(sim_time/dt+1); k++){
		re[i]->products[j]->delay_val[k] = (double*)malloc(sizeof(double)*4);
	      }
	    }
	    eq->number[index] = NULL;
	    eq->operator[index] = 0;
	    eq->delay_number[index] = re[i]->products[j]->delay_val;
	    eq->delay_comp_size[index] = NULL;
	    eq->explicit_delay_eq[index] = NULL;
	    if(initAssign != NULL){
	      for(k=0; k<num_of_time_variant_targets; k++){
		if(strcmp(time_variant_target_id[k], name) == 0){
		  for(k=0; k<Model_getNumInitialAssignments(m); k++){
		    if(strcmp(InitialAssignment_getSymbol(initAssign[k]->origin), name) == 0){
		      eq->explicit_delay_eq[index] = initAssign[k]->eq;
		    }
		  }
		}
	      }
	    }
	    if(timeVarAssign != NULL){
	      for(k=0; k<timeVarAssign->num_of_time_variant_assignments; k++){
		if(strcmp(timeVarAssign->target_id[k], name) == 0){
		  eq->explicit_delay_eq[index] = timeVarAssign->eq[k];
		}
	      }
	    }
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
	    //create delay
	    if(re[i]->reactants[j]->delay_val == NULL){
	      re[i]->reactants[j]->delay_val = (double**)malloc(sizeof(double*)*(int)(sim_time/dt+1));
	      for(k=0; k<(int)(sim_time/dt+1); k++){
		re[i]->reactants[j]->delay_val[k] = (double*)malloc(sizeof(double)*4);
	      }
	    }
	    eq->number[index] = NULL;
	    eq->operator[index] = 0;
	    eq->delay_number[index] = re[i]->reactants[j]->delay_val;
	    eq->delay_comp_size[index] = NULL;
	    eq->explicit_delay_eq[index] = NULL;
	    if(initAssign != NULL){
	      for(k=0; k<num_of_time_variant_targets; k++){
		if(strcmp(time_variant_target_id[k], name) == 0){
		  for(k=0; k<Model_getNumInitialAssignments(m); k++){
		    if(strcmp(InitialAssignment_getSymbol(initAssign[k]->origin), name) == 0){
		      eq->explicit_delay_eq[index] = initAssign[k]->eq;
		    }
		  }
		}
	      }
	    }
	    if(timeVarAssign != NULL){
	      for(k=0; k<timeVarAssign->num_of_time_variant_assignments; k++){
		if(strcmp(timeVarAssign->target_id[k], name) == 0){
		  eq->explicit_delay_eq[index] = timeVarAssign->eq[k];
		}
	      }
	    }
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
    index = get_equation(m, eq, sp, param, comp, re, left, index, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, mem);
  }
  if((right=ASTNode_getRightChild(node)) != NULL){
    index = get_equation(m, eq, sp, param, comp, re, right, index, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, mem);
  }
  if(ASTNode_isOperator(node)
     || ASTNode_isFunction(node)
     || ASTNode_isBoolean(node)){
    op = ASTNode_getType(node);
    eq->number[index] = NULL;
    eq->operator[index] = op;
    eq->delay_number[index] = NULL;
    eq->delay_comp_size[index] = NULL;
    eq->explicit_delay_eq[index] = NULL;
    index++;
  }else if(ASTNode_getType(node) == AST_NAME){
    name = ASTNode_getName(node);
    flag = 1;
    for(i=0; i<Model_getNumSpecies(m); i++){
      if(strcmp(name, Species_getId(sp[i]->origin)) == 0){
	eq->number[index] = &sp[i]->temp_value;
	eq->operator[index] = 0;
	eq->delay_number[index] = NULL;
	eq->delay_comp_size[index] = NULL;
	eq->explicit_delay_eq[index] = NULL;
	index++;
	flag = 0;
	break;
      }
    }
    if(flag){
      for(i=0; i<Model_getNumParameters(m); i++){
	if(strcmp(name, Parameter_getId(param[i]->origin)) == 0){
	  eq->number[index] = &param[i]->temp_value;
	  eq->operator[index] = 0;
	  eq->delay_number[index] = NULL;
	  eq->delay_comp_size[index] = NULL;
	  eq->explicit_delay_eq[index] = NULL;
	  index++;
	  flag = 0;
	  break;
	}
      }
    }
    if(flag){
      for(i=0; i<Model_getNumCompartments(m); i++){
	if(strcmp(name, Compartment_getId(comp[i]->origin)) == 0){
	  eq->number[index] = &comp[i]->temp_value;
	  eq->operator[index] = 0;
	  eq->delay_number[index] = NULL;
	  eq->delay_comp_size[index] = NULL;
	  eq->explicit_delay_eq[index] = NULL;
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
	    eq->number[index] = &re[i]->products[j]->temp_value;
	    eq->operator[index] = 0;
	    eq->delay_number[index] = NULL;
	    eq->delay_comp_size[index] = NULL;
	    eq->explicit_delay_eq[index] = NULL;
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
	    eq->number[index] = &re[i]->reactants[j]->temp_value;
	    eq->operator[index] = 0;
	    eq->delay_number[index] = NULL;
	    eq->delay_comp_size[index] = NULL;
	    eq->explicit_delay_eq[index] = NULL;
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
	eq->number[index] = time;
	eq->operator[index] = 0;
	eq->delay_number[index] = NULL;
	eq->delay_comp_size[index] = NULL;
	eq->explicit_delay_eq[index] = NULL;
	index++;
      }
    }
  }else if(ASTNode_getType(node) == AST_NAME_TIME){
    eq->number[index] = time;
    eq->operator[index] = 0;
    eq->delay_number[index] = NULL;
    eq->delay_comp_size[index] = NULL;
    eq->explicit_delay_eq[index] = NULL;
    index++;
  }else if(ASTNode_getType(node) == AST_NAME_AVOGADRO){
    eq->number[index] = (double*)malloc(sizeof(double));
    mem->memory[mem->num_of_allocated_memory++] = eq->number[index];
    value = 6.02214179e23;
    *eq->number[index] = value;
    eq->operator[index] = 0;
    eq->delay_number[index] = NULL;
    eq->delay_comp_size[index] = NULL;
    eq->explicit_delay_eq[index] = NULL;
    index++;
  }else if(ASTNode_getType(node) == AST_CONSTANT_E
	   || ASTNode_getType(node) == AST_CONSTANT_PI){
    eq->number[index] = (double*)malloc(sizeof(double));
    mem->memory[mem->num_of_allocated_memory++] = eq->number[index];
    if(ASTNode_getType(node) == AST_CONSTANT_E){
      value = 2.718281828459045235360287471352;
    }else{
      value = M_PI;
    }
    *eq->number[index] = value;
    eq->operator[index] = 0;
    eq->delay_number[index] = NULL;
    eq->delay_comp_size[index] = NULL;
    eq->explicit_delay_eq[index] = NULL;
    index++;
  }else{
    eq->number[index] = (double*)malloc(sizeof(double));
    mem->memory[mem->num_of_allocated_memory++] = eq->number[index];
    if(ASTNode_getType(node) == AST_INTEGER){
      value = ASTNode_getInteger(node);
    }else{
      value = ASTNode_getReal(node);
    }
    *eq->number[index] = value;
    eq->operator[index] = 0;
    eq->delay_number[index] = NULL;
    eq->delay_comp_size[index] = NULL;
    eq->explicit_delay_eq[index] = NULL;
    index++;
  }
  
  return index;
}
