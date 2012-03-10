#include "sbml/common/common.h"
#include "sbml/SBMLReader.h"
#include "sbml/SBMLTypes.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include "header.h"

void myASTNode_create(myASTNode *myNode, ASTNode_t *node, myASTNode *copied_myAST[], int *num_of_copied_myAST){
  ASTNode_t *left, *right;
  myASTNode *myLeft, *myRight;
  if((left=ASTNode_getLeftChild(node)) != NULL){
    myLeft = (myASTNode*)malloc(sizeof(myASTNode));
    copied_myAST[*num_of_copied_myAST] = myLeft;
    (*num_of_copied_myAST)++;
    myLeft->origin = left;
    myLeft->parent = myNode;
    myLeft->left = NULL;
    myLeft->right = NULL;
    myNode->left = myLeft;
    myASTNode_create(myLeft, left, copied_myAST, num_of_copied_myAST);
  }
  if((right=ASTNode_getRightChild(node)) != NULL){
    myRight = (myASTNode*)malloc(sizeof(myASTNode));
    copied_myAST[*num_of_copied_myAST] = myRight;
    (*num_of_copied_myAST)++;
    myRight->origin = right;
    myRight->parent = myNode;
    myRight->left = NULL;
    myRight->right = NULL;
    myNode->right = myRight;
    myASTNode_create(myRight, right, copied_myAST, num_of_copied_myAST);
  }
}

void ASTNode_recreate(myASTNode *myNode, ASTNode_t *node){
  if(myNode->left != NULL){
    ASTNode_replaceChild(node, 0, myNode->left->origin);
    ASTNode_recreate(myNode->left, ASTNode_getLeftChild(node));
  }
  if(myNode->right != NULL){
    ASTNode_replaceChild(node, 1, myNode->right->origin);
    ASTNode_recreate(myNode->right, ASTNode_getRightChild(node));
  }
}

void myASTNode_free(myASTNode *copied_myAST[], int num_of_copied_myAST){
  int i;
  for(i=0; i<num_of_copied_myAST; i++){
    free(copied_myAST[i]);
  }
}

void check_myAST(myASTNode *myNode){
  int type;
  if(myNode == NULL){
    return;
  }
  if(myNode->left != NULL){
    check_myAST(myNode->left);
  }
  if(myNode->right != NULL){
    check_myAST(myNode->right);
  }
  type = ASTNode_getType(myNode->origin);
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
    printf("integer(%ld) ", ASTNode_getInteger(myNode->origin));
    break;
  case AST_REAL:
    printf("real(%lf) ", ASTNode_getReal(myNode->origin));
    break;
  case AST_REAL_E:
    printf("real_E ");
    break;
  case AST_RATIONAL:
    printf("rational ");
    break;
  case AST_NAME:
    printf("name(%s) ", ASTNode_getName(myNode->origin));
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
    printf("function(%s) ", ASTNode_getName(myNode->origin));
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
  if(myNode->parent == NULL){
    printf("\n\n");
  }
}
