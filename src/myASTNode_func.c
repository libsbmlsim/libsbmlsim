#include "libsbmlsim/libsbmlsim.h"

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
    dbg_printf("+ ");
    break;
  case AST_MINUS:
    dbg_printf("- ");
    break;
  case AST_TIMES:
    dbg_printf("* ");
    break;
  case AST_DIVIDE:
    dbg_printf("/ ");
    break;
  case AST_POWER:
    dbg_printf("pow ");
    break;
  case AST_INTEGER:
    dbg_printf("integer(%ld) ", ASTNode_getInteger(myNode->origin));
    break;
  case AST_REAL:
    dbg_printf("real(%lf) ", ASTNode_getReal(myNode->origin));
    break;
  case AST_REAL_E:
    dbg_printf("real_E ");
    break;
  case AST_RATIONAL:
    dbg_printf("rational ");
    break;
  case AST_NAME:
    dbg_printf("name(%s) ", ASTNode_getName(myNode->origin));
    break;
  case AST_NAME_AVOGADRO:
    dbg_printf("avogadro ");
    break;
  case AST_NAME_TIME:
    dbg_printf("time ");
    break;
  case AST_CONSTANT_E:
    dbg_printf("constant ");
    break;
  case AST_CONSTANT_FALSE:
    dbg_printf("constant_false ");
    break;
  case AST_CONSTANT_PI:
    dbg_printf("pi ");
    break;
  case AST_CONSTANT_TRUE:
    dbg_printf("constant_true ");
    break;
  case AST_LAMBDA:
    dbg_printf("lambda ");
    break;
  case AST_FUNCTION:
    dbg_printf("function(%s) ", ASTNode_getName(myNode->origin));
    break;
  case AST_FUNCTION_ABS:
    dbg_printf("abs ");
    break;
  case AST_FUNCTION_ARCCOS:
    dbg_printf("arccos ");
    break;
  case AST_FUNCTION_ARCCOSH:
    dbg_printf("arccosh ");
    break;
  case AST_FUNCTION_ARCCOT:
    dbg_printf("arccot ");
    break;
  case AST_FUNCTION_ARCCOTH:
    dbg_printf("arccoth ");
    break;
  case AST_FUNCTION_ARCCSC:
    dbg_printf("arccsc ");
    break;
  case AST_FUNCTION_ARCCSCH:
    dbg_printf("arccsch ");
    break;
  case AST_FUNCTION_ARCSEC:
    dbg_printf("arcsec ");
    break;
  case AST_FUNCTION_ARCSECH:
    dbg_printf("arcsech ");
    break;
  case AST_FUNCTION_ARCSIN:
    dbg_printf("arcsin ");
    break;
  case AST_FUNCTION_ARCSINH:
    dbg_printf("arcsinh ");
    break;
  case AST_FUNCTION_ARCTAN:
    dbg_printf("arctan ");
    break;
  case AST_FUNCTION_ARCTANH:
    dbg_printf("arctanh ");
    break;
  case AST_FUNCTION_CEILING:
    dbg_printf("ceil ");
    break;
  case AST_FUNCTION_COS:
    dbg_printf("cos ");
    break;
  case AST_FUNCTION_COSH:
    dbg_printf("cosh ");
    break;
  case AST_FUNCTION_COT:
    dbg_printf("cot ");
    break;
  case AST_FUNCTION_COTH:
    dbg_printf("coth ");
    break;
  case AST_FUNCTION_CSC:
    dbg_printf("csc ");
    break;
  case AST_FUNCTION_CSCH:
    dbg_printf("csch ");
    break;
  case AST_FUNCTION_DELAY:
    dbg_printf("delay ");
    break;
  case AST_FUNCTION_EXP:
    dbg_printf("exp ");
    break;
  case AST_FUNCTION_FACTORIAL:
    dbg_printf("! ");
    break;
  case AST_FUNCTION_FLOOR:
    dbg_printf("floor ");
    break;
  case AST_FUNCTION_LN:
    dbg_printf("ln ");
    break;
  case AST_FUNCTION_LOG:
    dbg_printf("log10 ");
    break;
  case AST_FUNCTION_PIECEWISE:
    dbg_printf("piecewise ");
    break;
  case AST_FUNCTION_POWER:
    dbg_printf("f_pow ");
    break;
  case AST_FUNCTION_ROOT:
    dbg_printf("sqrt ");
    break;
  case AST_FUNCTION_SEC:
    dbg_printf("sec ");
    break;
  case AST_FUNCTION_SECH:
    dbg_printf("sech ");
    break;
  case AST_FUNCTION_SIN:
    dbg_printf("sin ");
    break;
  case AST_FUNCTION_SINH:
    dbg_printf("sinh ");
    break;
  case AST_FUNCTION_TAN:
    dbg_printf("tan ");
    break;
  case AST_FUNCTION_TANH:
    dbg_printf("tanh ");
    break;
  case AST_LOGICAL_AND:
    dbg_printf("and ");
    break;
  case AST_LOGICAL_NOT:
    dbg_printf("not ");
    break;
  case AST_LOGICAL_OR:
    dbg_printf("or ");
    break;
  case AST_LOGICAL_XOR:
    dbg_printf("xor ");
    break;
  case AST_RELATIONAL_EQ:
    dbg_printf("eq ");
    break;
  case AST_RELATIONAL_GEQ:
    dbg_printf("geq ");
    break;
  case AST_RELATIONAL_GT:
    dbg_printf("gt ");
    break;
  case AST_RELATIONAL_LEQ:
    dbg_printf("leq ");
    break;
  case AST_RELATIONAL_LT:
    dbg_printf("lt ");
    break;
  case AST_RELATIONAL_NEQ:
    dbg_printf("neq ");
    break;
  case AST_UNKNOWN:
    dbg_printf("unknown ");
    break;
  }
  if(myNode->parent == NULL){
    dbg_printf("\n\n");
  }
}
