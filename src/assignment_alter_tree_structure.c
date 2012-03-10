#include "sbml/common/common.h"
#include "sbml/SBMLReader.h"
#include "sbml/SBMLTypes.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include "header.h"

void assignment_alter_tree_structure(ASTNode_t **node_p, char* comp_name, int sw){
  ASTNode_t *times_node, *divide_node, *comp_node;
  if(sw == 0){
    times_node = ASTNode_createWithType(AST_TIMES);
    comp_node = ASTNode_createWithType(AST_NAME);
    ASTNode_setName(comp_node, comp_name);
    ASTNode_addChild(times_node, *node_p);
    ASTNode_addChild(times_node, comp_node);
    *node_p = times_node;
  }else if(sw == 1){
    divide_node = ASTNode_createWithType(AST_DIVIDE);
    comp_node = ASTNode_createWithType(AST_NAME);
    ASTNode_setName(comp_node, comp_name);
    ASTNode_addChild(divide_node, *node_p);
    ASTNode_addChild(divide_node, comp_node);
    *node_p = divide_node;    
  }
}

