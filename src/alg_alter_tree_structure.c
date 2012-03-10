#include "sbml/common/common.h"
#include "sbml/SBMLReader.h"
#include "sbml/SBMLTypes.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include "header.h"

void alg_alter_tree_structure(ASTNode_t **node_p, ASTNode_t *parent, int child_order){
  ASTNode_t *node, *left, *right;
  ASTNode_t *times_node, *one_node;
  node = *node_p;
  if((left=ASTNode_getLeftChild(node)) != NULL){
    alg_alter_tree_structure(&left, node, 0);
  }
  if((right=ASTNode_getRightChild(node)) != NULL){
    alg_alter_tree_structure(&right, node, 1);
  }
  if(ASTNode_getType(node) == AST_NAME){
    times_node = ASTNode_createWithType(AST_TIMES);
    one_node = ASTNode_createWithType(AST_INTEGER);
    ASTNode_setInteger(one_node, 1);
    ASTNode_addChild(times_node, one_node);
    ASTNode_addChild(times_node, node);
    if(parent != NULL){
      ASTNode_replaceChild(parent, child_order, times_node);
    }else{
      *node_p = times_node;
    }
  }
  return;
}
