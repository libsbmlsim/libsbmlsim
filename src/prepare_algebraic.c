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

void _prepare_algebraic1(ASTNode_t *node, char *included_id_in_alg[], unsigned int *num_of_included_id_in_alg);

void _prepare_algebraic2(Model_t *m, myASTNode *myNode, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], double sim_time, double dt, double *time, myInitialAssignment *initAssign[], char* time_variant_target_id[], unsigned int num_of_time_variant_targets, timeVariantAssignments *timeVarAssign, myAlgebraicEquations *algEq, int alg_order, char *target_id, int variable_order, allocated_memory *mem);

void _prepare_algebraic3(Model_t *m, ASTNode_t *node, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], double sim_time, double dt, double *time, myInitialAssignment *initAssign[], char *time_variant_target_id[], unsigned int num_of_time_variant_targets, timeVariantAssignments *timeVarAssign, myAlgebraicEquations *algEq, int alg_order, allocated_memory *mem);

void _prepare_algebraic4(ASTNode_t *node, myAlgebraicEquations *algEq);

/* find included id(species, parameter, compartment) in algebraic rule */
void _prepare_algebraic1(ASTNode_t *node, char *included_id_in_alg[], unsigned int *num_of_included_id_in_alg){
  unsigned int i;
  ASTNode_t *left, *right;
  int flag;
  left = ASTNode_getLeftChild(node);
  right = ASTNode_getRightChild(node);
  if(left != NULL){
    _prepare_algebraic1(left, included_id_in_alg, num_of_included_id_in_alg);
  }
  if(right != NULL){
    _prepare_algebraic1(right, included_id_in_alg, num_of_included_id_in_alg);
  }
  flag = 1;
  if(ASTNode_getType(node) == AST_NAME){
    for(i=0; i<*num_of_included_id_in_alg; i++){
      if(strcmp(ASTNode_getName(node), included_id_in_alg[i]) == 0){
        flag = 0;
      }
    }
    if(flag){
      included_id_in_alg[(*num_of_included_id_in_alg)++] = (char*)ASTNode_getName(node);
    }
  }
}

/* find coefficient tree */
void _prepare_algebraic2(Model_t *m, myASTNode *myNode, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], double sim_time, double dt, double *time, myInitialAssignment *initAssign[], char* time_variant_target_id[], unsigned int num_of_time_variant_targets, timeVariantAssignments *timeVarAssign, myAlgebraicEquations *algEq, int alg_order, char *target_id, int variable_order, allocated_memory *mem){
  ASTNode_t *minus_node, *zero_node, *final_eq_node;
  myASTNode *eq_root_node;
  int minus_sign;

  if(myNode->left != NULL){
    _prepare_algebraic2(m, myNode->left, sp, param, comp, re, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, algEq, alg_order, target_id, variable_order, mem);
  }
  if(myNode->right != NULL){
    _prepare_algebraic2(m, myNode->right, sp, param, comp, re, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, algEq, alg_order, target_id, variable_order, mem);
  }
  if(ASTNode_getType(myNode->origin) == AST_NAME){
    if(strcmp(ASTNode_getName(myNode->origin), target_id) == 0){
      ASTNode_setType(myNode->origin, AST_INTEGER);
      ASTNode_setInteger(myNode->origin, 1);
      eq_root_node = myNode;
      minus_sign = 1;
      while(eq_root_node->parent != NULL){
        if(ASTNode_getType(eq_root_node->parent->origin) != AST_TIMES
            && ASTNode_getType(eq_root_node->parent->origin) != AST_DIVIDE){
          if(ASTNode_getType(eq_root_node->parent->origin) == AST_MINUS
              && eq_root_node->parent->right == eq_root_node){
            minus_sign *= -1;
          }
          if(eq_root_node->parent->parent != NULL){
            if(eq_root_node->parent->parent->left == eq_root_node->parent){
              eq_root_node->parent->parent->left = eq_root_node;
            }else{
              eq_root_node->parent->parent->right = eq_root_node;
            }
            eq_root_node->parent = eq_root_node->parent->parent;
          }else{
            eq_root_node->parent = NULL;
            break;
          } 
        }else{
          eq_root_node = eq_root_node->parent;
        }
      }
      final_eq_node = eq_root_node->origin;
      TRACE(("myASTNode is\n"));
      check_myAST(eq_root_node);
      ASTNode_recreate(eq_root_node, final_eq_node);
      if(minus_sign == -1){
        minus_node = ASTNode_createWithType(AST_MINUS);
        zero_node = ASTNode_createWithType(AST_INTEGER);
        ASTNode_setInteger(zero_node, 0);
        ASTNode_addChild(minus_node, zero_node);
        ASTNode_addChild(minus_node, final_eq_node);
        final_eq_node = minus_node;
      }
      if(algEq->num_of_algebraic_variables > 1){
        TRACE(("math AST of coefficient matrix[%d][%d] is\n", alg_order, variable_order));
        check_AST(eq_root_node->origin, NULL);
        algEq->coefficient_matrix[alg_order][variable_order]->math_length = get_equation(m, algEq->coefficient_matrix[alg_order][variable_order], sp, param, comp, re, final_eq_node, 0, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, mem);
      }else{
        TRACE(("math AST of coefficient is\n"));
        check_AST(eq_root_node->origin, NULL);
        algEq->coefficient->math_length = get_equation(m, algEq->coefficient, sp, param, comp, re, final_eq_node, 0, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, mem);
      }
    }
  }
}

/* find constant vector */
void _prepare_algebraic3(Model_t *m, ASTNode_t *node, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], double sim_time, double dt, double *time, myInitialAssignment *initAssign[], char *time_variant_target_id[], unsigned int num_of_time_variant_targets, timeVariantAssignments *timeVarAssign, myAlgebraicEquations *algEq, int alg_order, allocated_memory *mem){
  _prepare_algebraic4(node, algEq);
  if(algEq->num_of_algebraic_variables > 1){
    TRACE(("math AST of constant vector[%d] is\n", alg_order));
    check_AST(node, NULL);
    algEq->constant_vector[alg_order]->math_length = get_equation(m, algEq->constant_vector[alg_order], sp, param, comp, re, node, 0, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, mem);
  }else{
    TRACE(("math AST of constant is\n"));
    check_AST(node, NULL);
    algEq->constant->math_length = get_equation(m, algEq->constant, sp, param, comp, re, node, 0, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, mem);
  }

}

/* recursive function in _prepare_algebraic3 */
void _prepare_algebraic4(ASTNode_t *node, myAlgebraicEquations *algEq){
  unsigned int i;
  ASTNode_t *left, *right;
  int flag;

  if((left=ASTNode_getLeftChild(node)) != NULL){
    _prepare_algebraic4(left, algEq);
  }
  if((right=ASTNode_getRightChild(node)) != NULL){
    _prepare_algebraic4(right, algEq);
  }
  if(ASTNode_getType(node) == AST_NAME){
    flag = 0;
    for(i=0; i<algEq->num_of_algebraic_variables; i++){
      if(strcmp(ASTNode_getName(node), algEq->variables_id[i]) == 0){
        flag = 1;
      }
    }
    if(flag){
      ASTNode_setType(node, AST_INTEGER);
      ASTNode_setInteger(node, 0);
    }
  }
  return;
}

void prepare_algebraic(Model_t *m, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], myRule *ru[], myEvent *ev[], myInitialAssignment *initAssign[], myAlgebraicEquations *algEq, double sim_time, double dt, double *time, char *time_variant_target_id[], unsigned int num_of_time_variant_targets, timeVariantAssignments *timeVarAssign, allocated_memory *mem, copied_AST *cp_AST){
  unsigned int i, j, k;
  char *constants_in_alg[MAX_ALGEBRAIC_CONSTANTS];
  char *included_id_in_alg[MAX_ALGEBRAIC_CONSTANTS];
  unsigned int num_of_constants_in_alg = 0;
  unsigned int num_of_included_id_in_alg = 0;
  int flag;
  Species_t *local_sp;
  Parameter_t *local_param;
  Compartment_t *local_comp;
  ASTNode_t *node;
  myASTNode *myNode = NULL;
  myASTNode *copied_myAST[MAX_COPIED_AST];
  unsigned int num_of_copied_myAST = 0;
  /* find constant in calculation algebraic rule */
  /* reaction target(reactants and products) */
  TRACE(("Reaction\n"));
  for(i=0; i<Model_getNumSpecies(m); i++){
    flag = 0;
    local_sp = (Species_t*)ListOf_get(Model_getListOfSpecies(m), i);
    for(j=0; j<Model_getNumReactions(m); j++){
      for(k=0; k<re[j]->num_of_products; k++){
        if(strcmp(Species_getId(re[j]->products[k]->mySp->origin), Species_getId(local_sp)) == 0){
          flag = 1;
        }
      }
      for(k=0; k<re[j]->num_of_reactants; k++){
        if(strcmp(Species_getId(re[j]->reactants[k]->mySp->origin), Species_getId(local_sp)) == 0){
          flag = 1;
        }
      }
    }
    if(flag){
      constants_in_alg[num_of_constants_in_alg++] = (char*)Species_getId(local_sp);
    }
  }
  /* rule target */
  TRACE(("Rule\n"));
  for(i=0; i<Model_getNumRules(m); i++){
    if(Rule_isRate(ru[i]->origin) || Rule_isAssignment(ru[i]->origin)){
      if(ru[i]->target_species != NULL){
        constants_in_alg[num_of_constants_in_alg++] = (char*)Species_getId(ru[i]->target_species->origin);
      }
      if(ru[i]->target_parameter != NULL){
        constants_in_alg[num_of_constants_in_alg++] = (char*)Parameter_getId(ru[i]->target_parameter->origin);
      }
      if(ru[i]->target_compartment != NULL){
        constants_in_alg[num_of_constants_in_alg++] = (char*)Compartment_getId(ru[i]->target_compartment->origin);
      }
    }
  }
  /* event target */
  TRACE(("Event\n"));
  for(i=0; i<Model_getNumEvents(m); i++){
    for(j=0; j<Event_getNumEventAssignments(ev[i]->origin); j++){
      flag = 1;
      for(k=0; k<num_of_constants_in_alg; k++){
        if(ev[i]->assignments[j]->target_species != NULL){
          if(strcmp(constants_in_alg[k], Species_getId(ev[i]->assignments[j]->target_species->origin)) == 0){
            flag = 0;
          }
        }
        if(ev[i]->assignments[j]->target_parameter != NULL){
          if(strcmp(constants_in_alg[k], Parameter_getId(ev[i]->assignments[j]->target_parameter->origin)) == 0){
            flag = 0;
          }	
        }
        if(ev[i]->assignments[j]->target_compartment != NULL){
          if(strcmp(constants_in_alg[k], Compartment_getId(ev[i]->assignments[j]->target_compartment->origin)) == 0){
            flag = 0;
          }
        }
      }
      if(flag){
        if(ev[i]->assignments[j]->target_species != NULL){
          constants_in_alg[num_of_constants_in_alg++] = (char*)Species_getId(ev[i]->assignments[j]->target_species->origin);
        }
        if(ev[i]->assignments[j]->target_species != NULL){
          constants_in_alg[num_of_constants_in_alg++] = (char*)Species_getId(ev[i]->assignments[j]->target_species->origin);
        }
        if(ev[i]->assignments[j]->target_species != NULL){
          constants_in_alg[num_of_constants_in_alg++] = (char*)Species_getId(ev[i]->assignments[j]->target_species->origin);
        }
      }
    }
  }
  /* initial assignment */
  TRACE(("Initial Assignment\n"));
  for(i=0; i<Model_getNumInitialAssignments(m); i++){
    flag = 1;
    for(j=0; j<num_of_constants_in_alg; j++){
      if(initAssign[i]->target_species != NULL){
        if(strcmp(constants_in_alg[i], Species_getId(initAssign[i]->target_species->origin)) == 0){
          flag = 0;
        }
      }
      if(initAssign[i]->target_parameter != NULL){
        if(strcmp(constants_in_alg[i], Parameter_getId(initAssign[i]->target_parameter->origin)) == 0){
          flag = 0;
        }	
      }
      if(initAssign[i]->target_compartment != NULL){
        if(strcmp(constants_in_alg[i], Compartment_getId(initAssign[i]->target_compartment->origin)) == 0){
          flag = 0;
        }
      }
    }
    if(flag){
      if(initAssign[i]->target_species != NULL){
        constants_in_alg[num_of_constants_in_alg++] = (char*)Species_getId(initAssign[i]->target_species->origin);
      }
      if(initAssign[i]->target_parameter != NULL){
        constants_in_alg[num_of_constants_in_alg++] = (char*)Parameter_getId(initAssign[i]->target_parameter->origin);
      }
      if(initAssign[i]->target_compartment != NULL){
        constants_in_alg[num_of_constants_in_alg++] = (char*)Compartment_getId(initAssign[i]->target_compartment->origin);
      }
    }
  }
  /* constant */
  TRACE(("Constant\n"));
  for(i=0; i<Model_getNumSpecies(m); i++){
    local_sp = (Species_t*)ListOf_get(Model_getListOfSpecies(m), i);
    if(Species_getConstant(local_sp)){
      constants_in_alg[num_of_constants_in_alg++] = (char*)Species_getId(local_sp);
    }
  }
  for(i=0; i<Model_getNumParameters(m); i++){
    local_param = (Parameter_t*)ListOf_get(Model_getListOfParameters(m), i);
    if(Parameter_getConstant(local_param)){
      constants_in_alg[num_of_constants_in_alg++] = (char*)Parameter_getId(local_param);
    }
  }
  for(i=0; i<Model_getNumCompartments(m); i++){
    local_comp = (Compartment_t*)ListOf_get(Model_getListOfCompartments(m), i);
    if(Compartment_getConstant(local_comp)){
      constants_in_alg[num_of_constants_in_alg++] = (char*)Compartment_getId(local_comp);
    }
  }

  TRACE(("constants in alg are\n"));
  for(i=0; i<num_of_constants_in_alg; i++){
    TRACE(("%s\n", constants_in_alg[i]));
  }

  for(i=0; i<Model_getNumRules(m); i++){
    if(ru[i]->is_algebraic){
      node = (ASTNode_t*)Rule_getMath(ru[i]->origin);
      node = ASTNode_deepCopy(node);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      _prepare_algebraic1(node, included_id_in_alg, &num_of_included_id_in_alg);
    }
  }

  TRACE(("included id in algebraic rules are\n"));
  for(i=0; i<num_of_included_id_in_alg; i++){
    TRACE(("%s\n", included_id_in_alg[i]));
  }

  for(i=0; i<num_of_included_id_in_alg; i++){
    flag = 1;
    for(j=0; j<num_of_constants_in_alg; j++){
      if(strcmp(included_id_in_alg[i], constants_in_alg[j]) == 0){
        flag = 0;
      }
    }
    if(flag){
      algEq->variables_id[algEq->num_of_algebraic_variables++] = (char*)included_id_in_alg[i];
    }
  }
  TRACE(("algebraic variable is\n"));
  for(i=0; i<algEq->num_of_algebraic_variables; i++){
    TRACE(("%s\n", algEq->variables_id[i]));
  }

  /* get coeffient */
  TRACE(("get coefficient matrix\n"));
  for(i=0; i<Model_getNumRules(m); i++){
    if(ru[i]->is_algebraic){
      node = (ASTNode_t*)Rule_getMath(ru[i]->origin);
      node = ASTNode_deepCopy(node);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      alg_alter_tree_structure(&node, NULL, 0);
      TRACE(("algebraic AST is\n"));
      check_AST(node, NULL);
      for(j=0; j<algEq->num_of_algebraic_variables; j++){
        myNode = (myASTNode*)malloc(sizeof(myASTNode));
        copied_myAST[num_of_copied_myAST++] = myNode;
        myNode->origin = node;
        myNode->parent = NULL;
        myNode->left = NULL;
        myNode->right = NULL;
        myASTNode_create(myNode, node, copied_myAST, &num_of_copied_myAST);
        _prepare_algebraic2(m, myNode, sp, param, comp, re, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, algEq, i, algEq->variables_id[j], j, mem);
        myASTNode_free(copied_myAST, num_of_copied_myAST);
      }
    }
  }
  TRACE(("zero substitute\n"));
  if(algEq->num_of_algebraic_rules > 1){
    for(i=0; i<algEq->num_of_algebraic_rules; i++){
      for(j=0; j<algEq->num_of_algebraic_rules; j++){
        if(algEq->coefficient_matrix[i][j]->math_length == 0){
          node = ASTNode_createWithType(AST_INTEGER);
          ASTNode_setInteger(node, 0);
          algEq->coefficient_matrix[i][j]->math_length = get_equation(m, algEq->coefficient_matrix[i][j], sp, param, comp, re, node, 0, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, mem);
        }
      }
    }
  }
  /* get constant */
  TRACE(("get constant\n"));
  for(i=0; i<Model_getNumRules(m); i++){
    if(ru[i]->is_algebraic){
      node = (ASTNode_t*)Rule_getMath(ru[i]->origin);
      node = ASTNode_deepCopy(node);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      _prepare_algebraic3(m, node, sp, param, comp, re, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, algEq, i, mem);
    }
  }

  return;
}
