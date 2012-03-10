#include "sbml/common/common.h"
#include "sbml/SBMLReader.h"
#include "sbml/SBMLTypes.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include <float.h>
#include "header.h"

#define EPSIRON 1.0e-8
#define DOUBLE_EQ(x, v) (((v-EPSIRON) < x) && (x < (v+EPSIRON)))

double calc(equation *eq, double dt, int cycle, double *reverse_time, int rk_order){
  int i;
  int pos = 0;
  int dummy = 0;
  double stack[eq->math_length];
  double **delay_preserver = NULL;
  double **delay_comp_preserver = NULL;
  double *delay_value = NULL;
  double *delay_comp_size = NULL;
  equation *explicit_delay_eq_preserver = NULL;

  for(i=0; i<eq->math_length; i++){
    if(eq->number[i]!=NULL){
      stack[pos] = *eq->number[i];
      //printf("%lf is stacked\n", stack[pos]);
      pos++;
    }else if(eq->delay_number[i]!=NULL){
      delay_preserver = eq->delay_number[i];
      delay_comp_preserver = eq->delay_comp_size[i];
      stack[pos] = dummy;
      if(eq->explicit_delay_eq[i]!=NULL){
	explicit_delay_eq_preserver = eq->explicit_delay_eq[i];
      }
      pos++;
    }else{
      switch(eq->operator[i]){
      case AST_PLUS: 
	//printf("operate +\n");
	stack[pos-2] += stack[pos-1];
	pos--;
	break;
      case AST_MINUS:
	//printf("operate -\n");
	stack[pos-2] -= stack[pos-1];
	pos--;
	break;
      case AST_TIMES:
	//printf("operate *\n");
	stack[pos-2] *= stack[pos-1];
	pos--;
	break;
      case AST_DIVIDE:
	//printf("operate /\n");
	stack[pos-2] /= stack[pos-1];
	pos--;
	break;
      case AST_POWER:
	//printf("operate pow\n");
	stack[pos-2] = pow(stack[pos-2], stack[pos-1]);
	pos--;
	break;
      case AST_FUNCTION_POWER:
	//printf("operate pow\n");
	stack[pos-2] = pow(stack[pos-2], stack[pos-1]);
	pos--;
	break;
      case AST_FUNCTION_FACTORIAL:
	//printf("operate !\n");
	stack[pos-1] = (double)factorial((int)stack[pos-1]);
	break;
      case AST_FUNCTION_ABS:	
	//printf("operate abs\n");
	stack[pos-1] = fabs(stack[pos-1]);
	break;
      case AST_FUNCTION_SIN:
	//printf("operate sin\n");
	stack[pos-1] = sin(stack[pos-1]);
	break;
      case AST_FUNCTION_COS:
	//printf("operate cos\n");
	stack[pos-1] = cos(stack[pos-1]);
	break;
      case AST_FUNCTION_TAN:
	//printf("operate tan\n");
	stack[pos-1] = tan(stack[pos-1]);
	break;
      case AST_FUNCTION_CSC:
	//printf("operate csc\n");
	stack[pos-1] = 1.0/sin(stack[pos-1]);
	break;
      case AST_FUNCTION_SEC:
	//printf("operate sec\n");
	stack[pos-1] = 1.0/cos(stack[pos-1]);
	break;
      case AST_FUNCTION_COT:
	//printf("operate cot\n");
	stack[pos-1] = 1.0/tan(stack[pos-1]);
	break;
      case AST_FUNCTION_ARCSIN:
	//printf("operate arcsin\n");
	if(stack[pos-1] > 1){
	  stack[pos-1] = asin(1);
	}else if(stack[pos-1] < -1){
	  stack[pos-1] = asin(-1);
	}else{
	  stack[pos-1] = asin(stack[pos-1]);
	}
	break;
      case AST_FUNCTION_ARCCOS:
	//printf("operate arccos\n");
	if(stack[pos-1] > 1){
	  stack[pos-1] = acos(1);
	}else if(stack[pos-1] < -1){
	  stack[pos-1] = acos(-1);
	}else{
	  stack[pos-1] = acos(stack[pos-1]);
	}
	break;
      case AST_FUNCTION_ARCTAN:
	//printf("operate arctan\n");
	stack[pos-1] = atan(stack[pos-1]);
	break;
      case AST_FUNCTION_ARCCSC:
	//printf("operate arccsc\n");
	if(1.0/stack[pos-1] > 1){
	  stack[pos-1] = asin(1);
	}else if(1.0/stack[pos-1] < -1){
	  stack[pos-1] = asin(-1);
	}else{
	  stack[pos-1] = asin(1.0/stack[pos-1]);
	}
	break;
      case AST_FUNCTION_ARCSEC:
	//printf("operate arcsec\n");
	if(1.0/stack[pos-1] > 1){
	  stack[pos-1] = acos(1);
	}else if(1.0/stack[pos-1] < -1){
	  stack[pos-1] = acos(-1);
	}else{
	  stack[pos-1] = acos(1.0/stack[pos-1]);
	}
	break;
      case AST_FUNCTION_ARCCOT:
	//printf("operate arccot\n");
	if(eq->operator[i-1] == AST_MINUS
	   && DOUBLE_EQ(*eq->number[i-2], 0)
	   && DOUBLE_EQ(*eq->number[i-3], 0)){
	  stack[pos-1] = atan(-1.0/stack[pos-1]);
	}else{
	  stack[pos-1] = atan(1.0/stack[pos-1]);
	}
	break;
      case AST_FUNCTION_SINH:
	//printf("operate sinh\n");
	stack[pos-1] = sinh(stack[pos-1]);
	break;
      case AST_FUNCTION_COSH:
	//printf("operate cosh\n");
	stack[pos-1] = cosh(stack[pos-1]);
	break;
      case AST_FUNCTION_TANH:
	//printf("operate tanh\n");
	stack[pos-1] = tanh(stack[pos-1]);
	break;
      case AST_FUNCTION_CSCH:
	//printf("operate csch\n");
	stack[pos-1] = sinh(1.0/stack[pos-1]);
	break;
      case AST_FUNCTION_SECH:
	//printf("operate sech\n");
	stack[pos-1] = cosh(1.0/stack[pos-1]);
	break;
      case AST_FUNCTION_COTH:
	//printf("operate coth\n");
	stack[pos-1] = tanh(1.0/stack[pos-1]);
	break;
      case AST_FUNCTION_ARCSINH:
	//printf("operate arcsinh\n");
	stack[pos-1] = asinh(stack[pos-1]);
	break;
      case AST_FUNCTION_ARCCOSH:
	//printf("operate arccosh\n");
	stack[pos-1] = acosh(stack[pos-1]);
	break;
      case AST_FUNCTION_ARCTANH:
	//printf("operate arctanh\n");
	if(stack[pos-1] >= 1){
	  stack[pos-1] = DBL_MAX;
	}else if(stack[pos-1] <= -1){
	  stack[pos-1] = -DBL_MAX;
	}else{
	  stack[pos-1] = atanh(stack[pos-1]);
	}
	break;
      case AST_FUNCTION_ARCCSCH:
	//printf("operate arccsch\n");
	stack[pos-1] = asinh(1.0/stack[pos-1]);
	break;
      case AST_FUNCTION_ARCSECH:
	//printf("operate arcsech\n");
	if(DOUBLE_EQ(stack[pos-1], 0)){
	  stack[pos-1] = DBL_MAX;
	}else if(stack[pos-1] > 1){
	  stack[pos-1] = 0;
	}else{
	  stack[pos-1] = acosh(1.0/stack[pos-1]);
	}
	break;
      case AST_FUNCTION_ARCCOTH:
	//printf("operate arccoth\n");
	if(1.0/stack[pos-1] >= 1){
	  stack[pos-1] = DBL_MAX;
	}else if(1.0/stack[pos-1] <= -1){
	  stack[pos-1] = -DBL_MAX;
	}else{
	  stack[pos-1] = atanh(1.0/stack[pos-1]);
	}
	break;
      case AST_FUNCTION_EXP:
	//printf("operate exp\n");
	stack[pos-1] = exp(stack[pos-1]);
	break;
      case AST_FUNCTION_LN:
	//printf("operate ln\n");
	stack[pos-1] = log(stack[pos-1]);
	break;
      case AST_FUNCTION_LOG:
	//printf("operate log\n");
	stack[pos-2] = log(stack[pos-1])/log(stack[pos-2]);
	pos--;
	break;
      case AST_FUNCTION_ROOT:
	//printf("operate root\n");
	stack[pos-2] = pow(stack[pos-1], 1/stack[pos-2]);
	pos--;
	break;
      case AST_FUNCTION_CEILING:
	//printf("operate ceiling\n");
	stack[pos-1] = ceilf((float)stack[pos-1]);
	break;
      case AST_FUNCTION_FLOOR:
	//printf("operate floor\n");
	stack[pos-1] = floorf((float)stack[pos-1]);
	break;
      case AST_FUNCTION_DELAY:
	//printf("operate delay\n");
	if(delay_comp_preserver != NULL){
	  if(cycle-(int)(stack[pos-1]/dt) > 0){
	    delay_value = delay_preserver[cycle-(int)(stack[pos-1]/dt)];
	    delay_comp_size = delay_comp_preserver[cycle-(int)(stack[pos-1]/dt)];
	    stack[pos-2] = delay_value[rk_order]/delay_comp_size[rk_order];
	  }else if(explicit_delay_eq_preserver != NULL){
	    *reverse_time = cycle*dt - stack[pos-1];
	    delay_comp_size = delay_comp_preserver[0];
       	    stack[pos-2] = calc(explicit_delay_eq_preserver, dt, cycle, reverse_time, rk_order)/delay_comp_size[rk_order];
	    explicit_delay_eq_preserver = NULL;
	  }else{
	    delay_value = delay_preserver[0];
	    delay_comp_size = delay_comp_preserver[0];
	    stack[pos-2] = delay_value[rk_order]/delay_comp_size[rk_order];
	  }
	}else{
	  if(cycle-(int)(stack[pos-1]/dt) > 0){
	    delay_value = delay_preserver[cycle-(int)(stack[pos-1]/dt)];
	    stack[pos-2] = delay_value[rk_order];
	  }else if(explicit_delay_eq_preserver != NULL){
	    *reverse_time = cycle*dt - stack[pos-1];
       	    stack[pos-2] = calc(explicit_delay_eq_preserver, dt, cycle, reverse_time, rk_order);
	    explicit_delay_eq_preserver = NULL;
	  }else{
	    delay_value = delay_preserver[0];
	    stack[pos-2] = delay_value[rk_order];
	  }
	}
	delay_preserver = NULL;
	delay_comp_preserver = NULL;
	pos--;
	break;
      case AST_RELATIONAL_EQ:
	//printf("operate eq\n");
	//if(stack[pos-2] == stack[pos-1]){
	//printf("EQ stack[pos-2] = %lf : stack[pos-1] = %lf\n", stack[pos-2], stack[pos-1]);
	if(DOUBLE_EQ(stack[pos-2], stack[pos-1])){
	  stack[pos-2] = 1;
	}else{
	  stack[pos-2] = 0;
	}
	pos--;
	break;
      case AST_RELATIONAL_NEQ:
	//printf("operate neq\n");
	//if(stack[pos-2] != stack[pos-1]){
	//printf("NEQ stack[pos-2] = %lf : stack[pos-1] = %lf\n", stack[pos-2], stack[pos-1]);
	if(!DOUBLE_EQ(stack[pos-2], stack[pos-1])){
	  stack[pos-2] = 1;
	}else{
	  stack[pos-2] = 0;
	}
	pos--;
	break;
      case AST_RELATIONAL_LT:
	//printf("operate lt\n");
	if(stack[pos-2] < stack[pos-1]){
	  stack[pos-2] = 1;
	}else{
	  stack[pos-2] = 0;
	}
	pos--;
	break;
      case AST_RELATIONAL_GT:
	//printf("operate gt\n");
	if(stack[pos-2] > stack[pos-1]){
	  stack[pos-2] = 1;
	}else{
	  stack[pos-2] = 0;
	}
	pos--;
	break;
      case AST_RELATIONAL_LEQ:
	//printf("operate leq\n");
	if(stack[pos-2] <= stack[pos-1]){
	  stack[pos-2] = 1;
	}else{
	  stack[pos-2] = 0;
	}
	pos--;
	break;
      case AST_RELATIONAL_GEQ:
	//printf("operate geq\n");
	if(stack[pos-2] >= stack[pos-1]){
	  stack[pos-2] = 1;
	}else{
	  stack[pos-2] = 0;
	}
	pos--;
	break;
      case AST_LOGICAL_AND:
	//printf("operate and\n");
	if(stack[pos-2] >= 0.5 && stack[pos-1] >= 0.5){
	  stack[pos-2] = 1;
	}else{
	  stack[pos-2] = 0;
	}
	pos--;
	break;
      case AST_LOGICAL_NOT:
	//printf("operate not\n");
	if(stack[pos-1] >= 0.5){
	  stack[pos-1] = 0;
	}else{
	  stack[pos-1] = 1;
	}
	break;
      case AST_LOGICAL_OR:
	//printf("operate or\n");
	if(stack[pos-2] >= 0.5 || stack[pos-1] >= 0.5){
	  stack[pos-2] = 1;
	}else{
	  stack[pos-2] = 0;
	}
	pos--;
	break;
      case AST_LOGICAL_XOR:
	//printf("operate xor\n");
	if((stack[pos-2] >= 0.5 && stack[pos-1] < 0.5)
	   || (stack[pos-2] < 0.5 && stack[pos-1] >= 0.5)){
	  stack[pos-2] = 1;
	}else{
	  stack[pos-2] = 0;
	}
	pos--;
	break;
      case AST_CONSTANT_TRUE:
	//printf("stack true\n");
	stack[pos] = 1;
	pos++;
	break;
      case AST_CONSTANT_FALSE:
	//printf("stack false\n");
	stack[pos] = 0;
	pos++;
	break;
      }
    }
  }

  return stack[0];
}
