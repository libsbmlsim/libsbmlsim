#include "sbml/common/common.h"
#include "sbml/SBMLReader.h"
#include "sbml/SBMLTypes.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include <unistd.h>

void check_AST(ASTNode_t* node, ASTNode_t* parent){
  int i;
  int type;
  ASTNode_t *checker;
  if(node == NULL){
    return;
  }
  checker = parent;
  for(i=0; i<ASTNode_getNumChildren(node); i++){
    parent = node;
    check_AST(ASTNode_getChild(node, i), parent);
  }
  type = ASTNode_getType(node);
  switch(type){
  case AST_PLUS:
    printf("+ ");
    break;
  case AST_MINUS:
    printf("- ");
    break;
  case AST_TIMES:
    printf("* ");
    break;
  case AST_DIVIDE:
    printf("/ ");
    break;
  case AST_POWER:
    printf("pow ");
    break;
  case AST_INTEGER:
    printf("integer(%ld) ", ASTNode_getInteger(node));
    break;
  case AST_REAL:
    printf("real(%lf) ", ASTNode_getReal(node));
    break;
  case AST_REAL_E:
    printf("real_E ");
    break;
  case AST_RATIONAL:
    printf("rational ");
    break;
  case AST_NAME:
    printf("name(%s) ", ASTNode_getName(node));
    break;
  case AST_NAME_AVOGADRO:
    printf("avogadro ");
    break;
  case AST_NAME_TIME:
    printf("time ");
    break;
  case AST_CONSTANT_E:
    printf("constant ");
    break;
  case AST_CONSTANT_FALSE:
    printf("constant_false ");
    break;
  case AST_CONSTANT_PI:
    printf("pi ");
    break;
  case AST_CONSTANT_TRUE:
    printf("constant_true ");
    break;
  case AST_LAMBDA:
    printf("lambda ");
    break;
  case AST_FUNCTION:
    printf("function(%s) ", ASTNode_getName(node));
    break;
  case AST_FUNCTION_ABS:
    printf("abs ");
    break;
  case AST_FUNCTION_ARCCOS:
    printf("arccos ");
    break;
  case AST_FUNCTION_ARCCOSH:
    printf("arccosh ");
    break;
  case AST_FUNCTION_ARCCOT:
    printf("arccot ");
    break;
  case AST_FUNCTION_ARCCOTH:
    printf("arccoth ");
    break;
  case AST_FUNCTION_ARCCSC:
    printf("arccsc ");
    break;
  case AST_FUNCTION_ARCCSCH:
    printf("arccsch ");
    break;
  case AST_FUNCTION_ARCSEC:
    printf("arcsec ");
    break;
  case AST_FUNCTION_ARCSECH:
    printf("arcsech ");
    break;
  case AST_FUNCTION_ARCSIN:
    printf("arcsin ");
    break;
  case AST_FUNCTION_ARCSINH:
    printf("arcsinh ");
    break;
  case AST_FUNCTION_ARCTAN:
    printf("arctan ");
    break;
  case AST_FUNCTION_ARCTANH:
    printf("arctanh ");
    break;
  case AST_FUNCTION_CEILING:
    printf("ceil ");
    break;
  case AST_FUNCTION_COS:
    printf("cos ");
    break;
  case AST_FUNCTION_COSH:
    printf("cosh ");
    break;
  case AST_FUNCTION_COT:
    printf("cot ");
    break;
  case AST_FUNCTION_COTH:
    printf("coth ");
    break;
  case AST_FUNCTION_CSC:
    printf("csc ");
    break;
  case AST_FUNCTION_CSCH:
    printf("csch ");
    break;
  case AST_FUNCTION_DELAY:
    printf("delay ");
    break;
  case AST_FUNCTION_EXP:
    printf("exp ");
    break;
  case AST_FUNCTION_FACTORIAL:
    printf("! ");
    break;
  case AST_FUNCTION_FLOOR:
    printf("floor ");
    break;
  case AST_FUNCTION_LN:
    printf("ln ");
    break;
  case AST_FUNCTION_LOG:
    printf("log10 ");
    break;
  case AST_FUNCTION_PIECEWISE:
    printf("piecewise ");
    break;
  case AST_FUNCTION_POWER:
    printf("f_pow ");
    break;
  case AST_FUNCTION_ROOT:
    printf("sqrt ");
    break;
  case AST_FUNCTION_SEC:
    printf("sec ");
    break;
  case AST_FUNCTION_SECH:
    printf("sech ");
    break;
  case AST_FUNCTION_SIN:
    printf("sin ");
    break;
  case AST_FUNCTION_SINH:
    printf("sinh ");
    break;
  case AST_FUNCTION_TAN:
    printf("tan ");
    break;
  case AST_FUNCTION_TANH:
    printf("tanh ");
    break;
  case AST_LOGICAL_AND:
    printf("and ");
    break;
  case AST_LOGICAL_NOT:
    printf("not ");
    break;
  case AST_LOGICAL_OR:
    printf("or ");
    break;
  case AST_LOGICAL_XOR:
    printf("xor ");
    break;
  case AST_RELATIONAL_EQ:
    printf("eq ");
    break;
  case AST_RELATIONAL_GEQ:
    printf("geq ");
    break;
  case AST_RELATIONAL_GT:
    printf("gt ");
    break;
  case AST_RELATIONAL_LEQ:
    printf("leq ");
    break;
  case AST_RELATIONAL_LT:
    printf("lt ");
    break;
  case AST_RELATIONAL_NEQ:
    printf("neq ");
    break;
  case AST_UNKNOWN:
    printf("unknown ");
    break;
  }
  if(checker == NULL){
    printf("\n\n");
  }
}
