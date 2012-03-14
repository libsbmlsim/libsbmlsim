#include "libsbmlsim/libsbmlsim.h"

void minus_func(ASTNode_t *node){
  int i;
  ASTNode_t *next_node, *zero_node;
  for(i=0; i<ASTNode_getNumChildren(node); i++){
    next_node = ASTNode_getChild(node, i);
    if(ASTNode_getNumChildren(node) == 1 && ASTNode_getType(node) == AST_MINUS){
      zero_node = ASTNode_create();
      ASTNode_setType(zero_node, AST_REAL);
      ASTNode_setReal(zero_node, 0);
      ASTNode_replaceChild(node, 0, zero_node);
      ASTNode_addChild(node, next_node);
    }else{
      minus_func(next_node);
    }
  }
  return;
}

void alter_tree_structure(Model_t *m, ASTNode_t **node_p, ASTNode_t *parent, int child_order, copied_AST *cp_AST){
  ASTNode_t *zero_node;
  ASTNode_t *compartment_node;
  ASTNode_t *node, *next_node;
  ASTNode_t *pc_eq, *pc_cd, *times_node, *and_node, *not_node, *divide_node;
  int i, j;
  ASTNode_t *arg_node_list[MAX_ARG_NUM];
  int arg_node_num;
  FunctionDefinition_t *fd;
  ASTNode_t *fd_arg;
  ASTNode_t *fd_body;
  Species_t *sp;

  node = *node_p;
  for(i=0; i<ASTNode_getNumChildren(node); i++){
    next_node = ASTNode_getChild(node, i);
    if(ASTNode_getNumChildren(node) == 1 && ASTNode_getType(node) == AST_MINUS){
      zero_node = ASTNode_create();
      ASTNode_setType(zero_node, AST_REAL);
      ASTNode_setReal(zero_node, 0);
      ASTNode_replaceChild(node, 0, zero_node);
      ASTNode_addChild(*node_p, next_node);
    }else{
      /* dbg_printf("down to %d th child from\n", i); */
      /* print_node_type(node); */
      alter_tree_structure(m, &next_node, *node_p, i, cp_AST);
    }
  }

  if(ASTNode_getType(node) == AST_NAME){
    for(i=0; i<Model_getNumSpecies(m); i++){
      sp = (Species_t*)ListOf_get(Model_getListOfSpecies(m), i);      
      if(strcmp(Species_getId(sp), ASTNode_getName(node)) == 0){
        if(!Species_getHasOnlySubstanceUnits(sp) && Species_isSetInitialAmount(sp) && Compartment_getSpatialDimensions(Model_getCompartmentById(m, Species_getCompartment(sp))) != 0){/* use val/comp in calculation */
          divide_node = ASTNode_createWithType(AST_DIVIDE);
          compartment_node = ASTNode_createWithType(AST_NAME);
          ASTNode_setName(compartment_node, Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(sp))));
          ASTNode_addChild(divide_node, node);
          ASTNode_addChild(divide_node, compartment_node);
          if(parent != NULL){
            ASTNode_replaceChild(parent, child_order, divide_node);
          }else{
            *node_p = divide_node;
          }
          node = *node_p;
          break;
        }else if(Species_getHasOnlySubstanceUnits(sp) && Species_isSetInitialConcentration(sp) && Compartment_getSpatialDimensions(Model_getCompartmentById(m, Species_getCompartment(sp))) != 0){/*  use val*comp in calculation */
          times_node = ASTNode_createWithType(AST_TIMES);
          compartment_node = ASTNode_createWithType(AST_NAME);
          ASTNode_setName(compartment_node, Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(sp))));
          ASTNode_addChild(times_node, node);
          ASTNode_addChild(times_node, compartment_node);
          if(parent != NULL){
            ASTNode_replaceChild(parent, child_order, times_node);
          }else{
            *node_p = times_node;
          }
          node = *node_p;
          break;
        }
      }
    }
  }
  if(ASTNode_getType(node) == AST_FUNCTION){
    arg_node_num = ASTNode_getNumChildren(node);
    for(i=0; i<arg_node_num; i++){
      arg_node_list[i] = ASTNode_getChild(node, i);
    }
    for(i=0; i<Model_getNumFunctionDefinitions(m); i++){
      fd = (FunctionDefinition_t*)ListOf_get(Model_getListOfFunctionDefinitions(m), i);
      fd_body = (ASTNode_t*)FunctionDefinition_getBody(fd);
      if(strcmp(FunctionDefinition_getId(fd), ASTNode_getName(node)) == 0){
        fd_body = ASTNode_deepCopy(fd_body);
        cp_AST->ast[cp_AST->num_of_copied_AST++] = fd_body;
        for(j=0; j<FunctionDefinition_getNumArguments(fd); j++){
          fd_arg = (ASTNode_t*)FunctionDefinition_getArgument(fd, j);
          ASTNode_replaceArgument(fd_body, (char*)ASTNode_getName(fd_arg), arg_node_list[j]);
        }
        /* test */
        /* 	for(i=0; i<ASTNode_getNumChildren(fd_body); i++){ */
        /* 	  next_node = ASTNode_getChild(fd_body, i); */
        /* 	  dbg_printf("down to %d th child from\n", i); */
        /* 	  print_node_type(node); */
        /* 	  alter_tree_structure(m, &next_node, fd_body, i, cp_AST); */
        /* 	} */
        minus_func(fd_body);
        check_AST(fd_body, NULL);
        if(parent != NULL){
          ASTNode_replaceChild(parent, child_order, fd_body);
        }else{
          *node_p = fd_body;
        }
        node = *node_p;
        break;
      }
    }
  }
  if(ASTNode_getType(node) == AST_FUNCTION_PIECEWISE){
    ASTNode_setType(node, AST_PLUS);
    ASTNode_setName(node, NULL);
    times_node = ASTNode_createWithType(AST_TIMES);
    pc_eq = ASTNode_getRightChild(node);
    ASTNode_addChild(times_node, pc_eq);
    if(ASTNode_getNumChildren(node) > 3){
      and_node = ASTNode_createWithType(AST_LOGICAL_AND);
      ASTNode_addChild(times_node, and_node);
      for(i=ASTNode_getNumChildren(node)-2; i >= 1; i = i-2){
        pc_cd = ASTNode_getChild(node, i);
        not_node = ASTNode_createWithType(AST_LOGICAL_NOT);
        ASTNode_addChild(not_node, pc_cd);
        ASTNode_addChild(and_node, not_node);
      }
      ASTNode_reduceToBinary(and_node);
    }else{
      pc_cd = ASTNode_getChild(node, 1);
      not_node = ASTNode_createWithType(AST_LOGICAL_NOT);
      ASTNode_addChild(not_node, pc_cd);
      ASTNode_addChild(times_node, not_node);
    }
    ASTNode_replaceChild(node, ASTNode_getNumChildren(node)-1, times_node);
    for(i=ASTNode_getNumChildren(node)-2; i >= 1; i = i-2){
      times_node = ASTNode_createWithType(AST_TIMES);
      pc_eq = ASTNode_getChild(node, i-1);
      pc_cd = ASTNode_getChild(node, i);
      ASTNode_addChild(times_node, pc_eq);
      ASTNode_addChild(times_node, pc_cd);
      ASTNode_removeChild(node, i);
      ASTNode_replaceChild(node ,i-1, times_node);
    }
    ASTNode_reduceToBinary(node);
  }
  /* print_node_type(node); */
  /* dbg_printf("is proccessed\n"); */
  return;
}
