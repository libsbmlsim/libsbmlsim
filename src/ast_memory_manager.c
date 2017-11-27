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
#include "libsbmlsim/libsbmlsim.h"

ast_memory_node_t* create_ast_memory_node(ASTNode_t* ast) {
#ifdef DEBUG_AST_MEMORY_DEBUG
  /* printf("create_ast_memory_node for [%p] %s\n", ast, SBML_formulaToString(ast)); */
#endif
  ast_memory_node_t* head = malloc(sizeof(ast_memory_node_t));
  head->ast = ast;
  head->next = NULL;
  return head;
}

void add_ast_memory_node(ASTNode_t* ast, char *file, int line) {
#ifdef DEBUG_AST_MEMORY_DEBUG
  /* printf("   add_ast_memory_node for (%s:%d) [%p] %s\n", file, line, ast, SBML_formulaToString(ast)); */
  printf("   add_ast_memory_node for (%s:%d) [%p]\n", file, line, ast);
#endif
  ast_memory_node_t* current;
  if (ast_memory_root_node == NULL) {
    ast_memory_root_node = create_ast_memory_node(ast);
  } else {
    current = ast_memory_root_node;
    while (current->next != NULL) {
      current = current->next;
    }
    current->next = create_ast_memory_node(ast);
  }
}

void remove_ast_memory_node(ASTNode_t* ast) {
#ifdef DEBUG_AST_MEMORY_DEBUG
  /* printf("remove_ast_memory_node for [%p] %s\n", ast, SBML_formulaToString(ast)); */
  printf("remove_ast_memory_node for [%p]\n", ast);
#endif
  ast_memory_node_t* tmp_node;
  ast_memory_node_t* current;
  if (ast_memory_root_node->ast == ast) {
    tmp_node = ast_memory_root_node;
    ast_memory_root_node = ast_memory_root_node->next;
    ASTNode_free(tmp_node->ast);
    free(tmp_node);
    return;
  }
  current = ast_memory_root_node;
  while (current->next != NULL) {
    if (current->next->ast == ast) {
      tmp_node = current->next;
      current->next = tmp_node->next;
      ASTNode_free(tmp_node->ast);
      free(tmp_node);
      return;
    }
    current = current->next;
  }
}

void free_all_ast_memory_nodes(void) {
  ast_memory_node_t *current = ast_memory_root_node;
  ast_memory_node_t *tmp_node;
  while (current != NULL) {
    tmp_node = current;
    current = tmp_node->next;
    remove_ast_memory_node(tmp_node->ast);
  }
}

