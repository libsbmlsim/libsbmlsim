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

void print_node_type(ASTNode_t *node){
  int type;
  if(node == NULL){
    return;
  }
  type = ASTNode_getType(node);
  switch(type){
    case AST_PLUS:
      printf("+\n");
      break;
    case AST_MINUS:
      printf("-\n");
      break;
    case AST_TIMES:
      printf("*\n");
      break;
    case AST_DIVIDE:
      printf("/\n");
      break;
    case AST_POWER:
      printf("pow\n");
      break;
    case AST_INTEGER:
      printf("integer(%ld)\n", ASTNode_getInteger(node));
      break;
    case AST_REAL:
      printf("real(%f)\n", ASTNode_getReal(node));
      break;
    case AST_REAL_E:
      printf("real_E\n");
      break;
    case AST_RATIONAL:
      printf("rational\n");
      break;
    case AST_NAME:
      printf("name(%s)\n", ASTNode_getName(node));
      break;
    case AST_NAME_AVOGADRO:
      printf("avogadro\n");
      break;
    case AST_NAME_TIME:
      printf("time\n");
      break;
    case AST_CONSTANT_E:
      printf("constant\n");
      break;
    case AST_CONSTANT_FALSE:
      printf("constant_false\n");
      break;
    case AST_CONSTANT_PI:
      printf("pi\n");
      break;
    case AST_CONSTANT_TRUE:
      printf("constant_true\n");
      break;
    case AST_LAMBDA:
      printf("lambda\n");
      break;
    case AST_FUNCTION:
      printf("function(%s)\n", ASTNode_getName(node));
      break;
    case AST_FUNCTION_ABS:
      printf("abs\n");
      break;
    case AST_FUNCTION_ARCCOS:
      printf("arccos\n");
      break;
    case AST_FUNCTION_ARCCOSH:
      printf("arccosh\n");
      break;
    case AST_FUNCTION_ARCCOT:
      printf("arccot\n");
      break;
    case AST_FUNCTION_ARCCOTH:
      printf("arccoth\n");
      break;
    case AST_FUNCTION_ARCCSC:
      printf("arccsc\n");
      break;
    case AST_FUNCTION_ARCCSCH:
      printf("arccsch\n");
      break;
    case AST_FUNCTION_ARCSEC:
      printf("arcsec\n");
      break;
    case AST_FUNCTION_ARCSECH:
      printf("arcsech\n");
      break;
    case AST_FUNCTION_ARCSIN:
      printf("arcsin\n");
      break;
    case AST_FUNCTION_ARCSINH:
      printf("arcsinh\n");
      break;
    case AST_FUNCTION_ARCTAN:
      printf("arctan\n");
      break;
    case AST_FUNCTION_ARCTANH:
      printf("arctanh\n");
      break;
    case AST_FUNCTION_CEILING:
      printf("ceil\n");
      break;
    case AST_FUNCTION_COS:
      printf("cos\n");
      break;
    case AST_FUNCTION_COSH:
      printf("cosh\n");
      break;
    case AST_FUNCTION_COT:
      printf("cot\n");
      break;
    case AST_FUNCTION_COTH:
      printf("coth\n");
      break;
    case AST_FUNCTION_CSC:
      printf("csc\n");
      break;
    case AST_FUNCTION_CSCH:
      printf("csch\n");
      break;
    case AST_FUNCTION_DELAY:
      printf("delay\n");
      break;
    case AST_FUNCTION_EXP:
      printf("exp\n");
      break;
    case AST_FUNCTION_FACTORIAL:
      printf("!\n");
      break;
    case AST_FUNCTION_FLOOR:
      printf("floor\n");
      break;
    case AST_FUNCTION_LN:
      printf("ln\n");
      break;
    case AST_FUNCTION_LOG:
      printf("log10\n");
      break;
    case AST_FUNCTION_PIECEWISE:
      printf("piecewise\n");
      break;
    case AST_FUNCTION_POWER:
      printf("f_pow\n");
      break;
    case AST_FUNCTION_ROOT:
      printf("sqrt\n");
      break;
    case AST_FUNCTION_SEC:
      printf("sec\n");
      break;
    case AST_FUNCTION_SECH:
      printf("sech\n");
      break;
    case AST_FUNCTION_SIN:
      printf("sin\n");
      break;
    case AST_FUNCTION_SINH:
      printf("sinh\n");
      break;
    case AST_FUNCTION_TAN:
      printf("tan\n");
      break;
    case AST_FUNCTION_TANH:
      printf("tanh\n");
      break;
    case AST_LOGICAL_AND:
      printf("and\n");
      break;
    case AST_LOGICAL_NOT:
      printf("not\n");
      break;
    case AST_LOGICAL_OR:
      printf("or\n");
      break;
    case AST_LOGICAL_XOR:
      printf("xor\n");
      break;
    case AST_RELATIONAL_EQ:
      printf("eq\n");
      break;
    case AST_RELATIONAL_GEQ:
      printf("geq\n");
      break;
    case AST_RELATIONAL_GT:
      printf("gt\n");
      break;
    case AST_RELATIONAL_LEQ:
      printf("leq\n");
      break;
    case AST_RELATIONAL_LT:
      printf("lt\n");
      break;
    case AST_RELATIONAL_NEQ:
      printf("neq\n");
      break;
    case AST_UNKNOWN:
      printf("unknown\n");
      break;
  }
}
