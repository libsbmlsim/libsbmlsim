#include "libsbmlsim/libsbmlsim.h"

void _prepare_reversible_fast_reaction(Model_t *m, myASTNode *myNode, myReaction *re, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re_whole[], double sim_time, double dt, double *time, myInitialAssignment *initAssign[], char *time_variant_target_id[], int num_of_time_variant_targets, timeVariantAssignments *timeVarAssign, char *target_id, int p_or_r, allocated_memory *mem){
  ASTNode_t *minus_node, *zero_node, *final_eq_node;
  myASTNode *eq_root_node;
  int minus_sign;
  
  if(myNode->left != NULL){
    _prepare_reversible_fast_reaction(m, myNode->left, re, sp, param, comp, re_whole, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, target_id, p_or_r, mem);
  }
  if(myNode->right != NULL){
    _prepare_reversible_fast_reaction(m, myNode->right, re, sp, param, comp, re_whole, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, target_id, p_or_r, mem);
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
      dbg_printf("myASTNode is\n");
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
      if(p_or_r == 0){//products coefficient
	dbg_printf("AST of product numerator is\n");
	check_AST(final_eq_node, NULL);
	re->products_equili_numerator->math_length = get_equation(m, re->products_equili_numerator, sp, param, comp, re_whole, final_eq_node, 0, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, mem);
      }else{//reactants coefficient
	minus_node = ASTNode_createWithType(AST_MINUS);
	zero_node = ASTNode_createWithType(AST_INTEGER);
	ASTNode_setInteger(zero_node, 0);
	ASTNode_addChild(minus_node, zero_node);
	ASTNode_addChild(minus_node, final_eq_node);
	final_eq_node = minus_node;
	dbg_printf("AST of reactant numerator is\n");
	check_AST(final_eq_node, NULL);
	re->reactants_equili_numerator->math_length = get_equation(m, re->reactants_equili_numerator, sp, param, comp, re_whole, final_eq_node, 0, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, mem);
      }
    }
  }
}

void prepare_reversible_fast_reaction(Model_t *m, myReaction *re[], mySpecies *sp[], myParameter *param[], myCompartment *comp[], double sim_time, double dt, double *time, myInitialAssignment *initAssign[], char *time_variant_target_id[], int num_of_time_variant_targets, timeVariantAssignments *timeVarAssign, allocated_memory *mem, copied_AST *cp_AST){
  int i;
  int num_of_reactions = Model_getNumReactions(m);
  ASTNode_t *node, *cp_node1, *cp_node2;
  myASTNode *myNode = NULL;
  myASTNode *copied_myAST[MAX_COPIED_AST];
  int num_of_copied_myAST = 0;
  for(i=0; i<num_of_reactions; i++){
    if(re[i]->is_fast && re[i]->is_reversible){
      node = (ASTNode_t*)KineticLaw_getMath(Reaction_getKineticLaw(re[i]->origin));
      node = ASTNode_deepCopy(node);
      dbg_printf("original math of %s: ", Reaction_getId(re[i]->origin));
      check_AST(node, NULL);
      //alter_tree_structure(m, &node, cp_AST);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      set_local_para_as_value(node, Reaction_getKineticLaw(re[i]->origin));
      dbg_printf("alterated math of %s : ", Reaction_getId(re[i]->origin));
      check_AST(node, NULL);
      cp_node1 = ASTNode_deepCopy(node);
      cp_node2 = ASTNode_deepCopy(node);
      //get products numerator
      myNode = (myASTNode*)malloc(sizeof(myASTNode));
      copied_myAST[num_of_copied_myAST++] = myNode;
      myNode->origin = cp_node1;
      myNode->parent = NULL;
      myNode->left = NULL;
      myNode->right = NULL;
      myASTNode_create(myNode, cp_node1, copied_myAST, &num_of_copied_myAST);
      re[i]->products_equili_numerator = (equation*)malloc(sizeof(equation));
      dbg_printf("target_id is %s\n", Species_getId(re[i]->reactants[0]->mySp->origin));
      check_AST(cp_node1, NULL);
      _prepare_reversible_fast_reaction(m, myNode, re[i], sp, param, comp, re, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, (char*)Species_getId(re[i]->reactants[0]->mySp->origin), 0, mem);
      //get reactants numerator
      myNode = (myASTNode*)malloc(sizeof(myASTNode));
      copied_myAST[num_of_copied_myAST++] = myNode;
      myNode->origin = cp_node2;
      myNode->parent = NULL;
      myNode->left = NULL;
      myNode->right = NULL;
      re[i]->reactants_equili_numerator = (equation*)malloc(sizeof(equation));
      myASTNode_create(myNode, cp_node2, copied_myAST, &num_of_copied_myAST);
      dbg_printf("target_id is %s\n", Species_getId(re[i]->products[0]->mySp->origin));
      check_AST(cp_node2, NULL);
      _prepare_reversible_fast_reaction(m, myNode, re[i], sp, param, comp, re, sim_time, dt, time, initAssign, time_variant_target_id, num_of_time_variant_targets, timeVarAssign, (char*)Species_getId(re[i]->products[0]->mySp->origin), 1, mem);
      myASTNode_free(copied_myAST, num_of_copied_myAST);
    }
  }
}
