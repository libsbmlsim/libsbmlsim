#include "libsbmlsim/libsbmlsim.h"

void set_local_para_as_value(ASTNode_t *node, KineticLaw_t *kl){
  int i;
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
