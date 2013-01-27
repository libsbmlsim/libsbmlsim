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
#include "libsbmlsim/copied_AST.h"
#include <stdlib.h>
#include <sbml/SBMLTypes.h>

copied_AST *copied_AST_create() {
  copied_AST *ast = (copied_AST *)malloc(sizeof(copied_AST));
  ast->num_of_copied_AST = 0;
  return ast;
}

void copied_AST_free(copied_AST *ast) {
  unsigned int i;

  if (ast == NULL) {
    return;
  }

  for (i = 0; i < ast->num_of_copied_AST; i++) {
    ASTNode_free(ast->ast[i]);
  }
  free(ast);
}

