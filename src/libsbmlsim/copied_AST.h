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
#ifndef LibSBMLSim_CopiedAST_h
#define LibSBMLSim_CopiedAST_h

#include "typedefs.h"
#include "common.h"
#include <sbml/SBMLTypes.h>

struct _copied_AST {
  ASTNode_t *ast[MAX_COPIED_AST];
  unsigned int num_of_copied_AST;
};

copied_AST *copied_AST_create();
void copied_AST_free(copied_AST *ast);

#endif /* LibSBMLSim_CopiedAST_h */
