/**
 * <!--------------------------------------------------------------------------
 * This file is part of libSBMLSim.  Please visit
 * http://fun.bio.keio.ac.jp/Software/libsbmlsim/ for more
 * information about libSBMLSim and its latest version.
 *
 * Copyright (C) 2011-2012 by the Keio University, Yokohama, Japan
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution.
 * ---------------------------------------------------------------------- -->*/
#include "libsbmlsim/libsbmlsim.h"

void set_local_para_as_value(ASTNode_t *node, KineticLaw_t *kl){
  unsigned int i;
  double value;
  ASTNode_t *left, *right;
  const char *name, *id;

  if((left=ASTNode_getLeftChild(node)) != NULL){
    set_local_para_as_value(left, kl);
  }
  if((right=ASTNode_getRightChild(node)) != NULL){
    set_local_para_as_value(right, kl);
  }
  if(ASTNode_getType(node) == AST_NAME){
    name = ASTNode_getName(node);
    for(i=0;i<KineticLaw_getNumParameters(kl);i++){
      id = Parameter_getId((Parameter_t*)ListOf_get(KineticLaw_getListOfParameters(kl), i));
      if(strcmp(name, id) == 0){
        value = Parameter_getValue((Parameter_t*)ListOf_get(KineticLaw_getListOfParameters(kl), i));	
        ASTNode_setType(node, AST_REAL);
        ASTNode_setReal(node, value);
      }
    }	
  }else if(ASTNode_getType(node) == AST_INTEGER){
    value = (double)ASTNode_getInteger(node);
    ASTNode_setType(node, AST_REAL);
    ASTNode_setReal(node, value);
  }
} 
