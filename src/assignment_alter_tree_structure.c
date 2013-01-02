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

