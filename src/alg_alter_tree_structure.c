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
