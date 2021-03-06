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

void check_math(equation *eq){
  unsigned int i;
  TRACE(("value:"));
  for(i=0; i<eq->math_length; i++){
    if(eq->number[i] != NULL){
      TRACE(("%lf ", *eq->number[i]));
    }else{
      TRACE(("NULL "));
    }
  }
  TRACE(("\n"));
  TRACE(("operator:"));
  for(i=0; i<eq->math_length; i++){
    switch(eq->op[i]){
      case AST_PLUS:
        TRACE(("+ "));
        break;
      case AST_MINUS:
        TRACE(("- "));
        break;
      case AST_TIMES:
        TRACE(("* "));
        break;
      case AST_DIVIDE:
        TRACE(("/ "));
        break;
      case AST_POWER:
        TRACE(("pow "));
        break;
      case AST_INTEGER:
        TRACE(("integer "));
        break;
      case AST_REAL:
        TRACE(("real "));
        break;
      case AST_REAL_E:
        TRACE(("real_E "));
        break;
      case AST_RATIONAL:
        TRACE(("rational "));
        break;
      case AST_NAME:
        TRACE(("name "));
        break;
      case AST_NAME_AVOGADRO:
        TRACE(("avogadro "));
        break;
      case AST_NAME_TIME:
        TRACE(("time "));
        break;
      case AST_CONSTANT_E:
        TRACE(("constant "));
        break;
      case AST_CONSTANT_FALSE:
        TRACE(("constant_false "));
        break;
      case AST_CONSTANT_PI:
        TRACE(("pi "));
        break;
      case AST_CONSTANT_TRUE:
        TRACE(("constant_true "));
        break;
      case AST_LAMBDA:
        TRACE(("lambda "));
        break;
      case AST_FUNCTION:
        TRACE(("function "));
        break;
      case AST_FUNCTION_ABS:
        TRACE(("abs "));
        break;
      case AST_FUNCTION_ARCCOS:
        TRACE(("arccos "));
        break;
      case AST_FUNCTION_ARCCOSH:
        TRACE(("arccosh "));
        break;
      case AST_FUNCTION_ARCCOT:
        TRACE(("arccot "));
        break;
      case AST_FUNCTION_ARCCOTH:
        TRACE(("arccoth "));
        break;
      case AST_FUNCTION_ARCCSC:
        TRACE(("arccsc "));
        break;
      case AST_FUNCTION_ARCCSCH:
        TRACE(("arccsch "));
        break;
      case AST_FUNCTION_ARCSEC:
        TRACE(("arcsec "));
        break;
      case AST_FUNCTION_ARCSECH:
        TRACE(("arcsech "));
        break;
      case AST_FUNCTION_ARCSIN:
        TRACE(("arcsin "));
        break;
      case AST_FUNCTION_ARCSINH:
        TRACE(("arcsinh "));
        break;
      case AST_FUNCTION_ARCTAN:
        TRACE(("arctan "));
        break;
      case AST_FUNCTION_ARCTANH:
        TRACE(("arctanh "));
        break;
      case AST_FUNCTION_CEILING:
        TRACE(("ceil "));
        break;
      case AST_FUNCTION_COS:
        TRACE(("cos "));
        break;
      case AST_FUNCTION_COSH:
        TRACE(("cosh "));
        break;
      case AST_FUNCTION_COT:
        TRACE(("cot "));
        break;
      case AST_FUNCTION_COTH:
        TRACE(("coth "));
        break;
      case AST_FUNCTION_CSC:
        TRACE(("csc "));
        break;
      case AST_FUNCTION_CSCH:
        TRACE(("csch "));
        break;
      case AST_FUNCTION_DELAY:
        TRACE(("delay "));
        break;
      case AST_FUNCTION_EXP:
        TRACE(("exp "));
        break;
      case AST_FUNCTION_FACTORIAL:
        TRACE(("! "));
        break;
      case AST_FUNCTION_FLOOR:
        TRACE(("floor "));
        break;
      case AST_FUNCTION_LN:
        TRACE(("ln "));
        break;
      case AST_FUNCTION_LOG:
        TRACE(("log10 "));
        break;
      case AST_FUNCTION_PIECEWISE:
        TRACE(("piecewise "));
        break;
      case AST_FUNCTION_POWER:
        TRACE(("f_pow "));
        break;
      case AST_FUNCTION_ROOT:
        TRACE(("sqrt "));
        break;
      case AST_FUNCTION_SEC:
        TRACE(("sec "));
        break;
      case AST_FUNCTION_SECH:
        TRACE(("sech "));
        break;
      case AST_FUNCTION_SIN:
        TRACE(("sin "));
        break;
      case AST_FUNCTION_SINH:
        TRACE(("sinh "));
        break;
      case AST_FUNCTION_TAN:
        TRACE(("tan "));
        break;
      case AST_FUNCTION_TANH:
        TRACE(("tanh "));
        break;
      case AST_LOGICAL_AND:
        TRACE(("and "));
        break;
      case AST_LOGICAL_NOT:
        TRACE(("not "));
        break;
      case AST_LOGICAL_OR:
        TRACE(("or "));
        break;
      case AST_LOGICAL_XOR:
        TRACE(("xor "));
        break;
      case AST_RELATIONAL_EQ:
        TRACE(("eq "));
        break;
      case AST_RELATIONAL_GEQ:
        TRACE(("geq "));
        break;
      case AST_RELATIONAL_GT:
        TRACE(("gt "));
        break;
      case AST_RELATIONAL_LEQ:
        TRACE(("leq "));
        break;
      case AST_RELATIONAL_LT:
        TRACE(("lt "));
        break;
      case AST_RELATIONAL_NEQ:
        TRACE(("neq "));
        break;
      case AST_UNKNOWN:
        TRACE(("unknown "));
        break;
      default:
        TRACE(("0 "));
    }
  }
  TRACE(("\n"));
  TRACE(("delay_value:"));
  for(i=0; i<eq->math_length; i++){
    if(eq->delay_number[i] != NULL){
      TRACE(("exist "));
    }else{
      TRACE(("NULL "));
    }
  }
  TRACE(("\n\n"));
}
