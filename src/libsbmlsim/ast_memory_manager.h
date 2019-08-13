/**
 * <!--------------------------------------------------------------------------
 * This file is part of libSBMLSim.  Please visit
 * http://fun.bio.keio.ac.jp/software/libsbmlsim/ for more
 * information about libSBMLSim and its latest version.
 *
 * Copyright (C) 2011-2017 by the Keio University, Yokohama, Japan
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution.
 * ---------------------------------------------------------------------- -->*/
#ifndef LibSBMLSim_ASTMemoryManager_h
#define LibSBMLSim_ASTMemoryManager_h

#include "common.h"

/* linked list */
typedef struct ast_memory_node {ASTNode_t* ast; struct ast_memory_node* next; } ast_memory_node_t;

ast_memory_node_t* create_ast_memory_node(ASTNode_t*);
void add_ast_memory_node(ASTNode_t*, char*, int);
void remove_ast_memory_node(ASTNode_t*);
void free_all_ast_memory_nodes(void);

#endif /* LibSBMLSim_ASTMemoryManager_h */
