#include "libsbmlsim/libsbmlsim.h"

double calc(equation *eq, double dt, int cycle, double *reverse_time, int rk_order){
  int i;
  int pos = 0;
  int dummy = 0;
  double **delay_preserver = NULL;
  double **delay_comp_preserver = NULL;
  double *delay_value = NULL;
  double *delay_comp_size = NULL;
  equation *explicit_delay_eq_preserver = NULL;
  /* double stack[eq->math_length]; */
  double *stack;

  /* XXX must free? */
  stack = (double*)malloc(sizeof(double) * eq->math_length);

  for(i=0; i<eq->math_length; i++){
    if(eq->number[i]!=NULL){
      stack[pos] = *eq->number[i];
      /* TRACE(("%lf is stacked\n", stack[pos])); */
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
          /* TRACE(("operate +\n")); */
          stack[pos-2] += stack[pos-1];
          pos--;
          break;
        case AST_MINUS:
          /* TRACE(("operate -\n")); */
          stack[pos-2] -= stack[pos-1];
          pos--;
          break;
        case AST_TIMES:
          /* TRACE(("operate *\n")); */
          stack[pos-2] *= stack[pos-1];
          pos--;
          break;
        case AST_DIVIDE:
          /* TRACE(("operate /\n")); */
          stack[pos-2] /= stack[pos-1];
          pos--;
          break;
        case AST_POWER:
          /* TRACE(("operate pow\n")); */
          stack[pos-2] = pow(stack[pos-2], stack[pos-1]);
          pos--;
          break;
        case AST_FUNCTION_POWER:
          /* TRACE(("operate pow\n")); */
          stack[pos-2] = pow(stack[pos-2], stack[pos-1]);
          pos--;
          break;
        case AST_FUNCTION_FACTORIAL:
          /* TRACE(("operate !\n")); */
          stack[pos-1] = (double)factorial((int)stack[pos-1]);
          break;
        case AST_FUNCTION_ABS:	
          /* TRACE(("operate abs\n")); */
          stack[pos-1] = fabs(stack[pos-1]);
          break;
        case AST_FUNCTION_SIN:
          /* TRACE(("operate sin\n")); */
          stack[pos-1] = sin(stack[pos-1]);
          break;
        case AST_FUNCTION_COS:
          /* TRACE(("operate cos\n")); */
          stack[pos-1] = cos(stack[pos-1]);
          break;
        case AST_FUNCTION_TAN:
          /* TRACE(("operate tan\n")); */
          stack[pos-1] = tan(stack[pos-1]);
          break;
        case AST_FUNCTION_CSC:
          /* TRACE(("operate csc\n")); */
          stack[pos-1] = 1.0/sin(stack[pos-1]);
          break;
        case AST_FUNCTION_SEC:
          /* TRACE(("operate sec\n")); */
          stack[pos-1] = 1.0/cos(stack[pos-1]);
          break;
        case AST_FUNCTION_COT:
          /* TRACE(("operate cot\n")); */
          stack[pos-1] = 1.0/tan(stack[pos-1]);
          break;
        case AST_FUNCTION_ARCSIN:
          /* TRACE(("operate arcsin\n")); */
          if(stack[pos-1] > 1){
            stack[pos-1] = asin(1);
          }else if(stack[pos-1] < -1){
            stack[pos-1] = asin(-1);
          }else{
            stack[pos-1] = asin(stack[pos-1]);
          }
          break;
        case AST_FUNCTION_ARCCOS:
          /* TRACE(("operate arccos\n")); */
          if(stack[pos-1] > 1){
            stack[pos-1] = acos(1);
          }else if(stack[pos-1] < -1){
            stack[pos-1] = acos(-1);
          }else{
            stack[pos-1] = acos(stack[pos-1]);
          }
          break;
        case AST_FUNCTION_ARCTAN:
          /* TRACE(("operate arctan\n")); */
          stack[pos-1] = atan(stack[pos-1]);
          break;
        case AST_FUNCTION_ARCCSC:
          /* TRACE(("operate arccsc\n")); */
          if(1.0/stack[pos-1] > 1){
            stack[pos-1] = asin(1);
          }else if(1.0/stack[pos-1] < -1){
            stack[pos-1] = asin(-1);
          }else{
            stack[pos-1] = asin(1.0/stack[pos-1]);
          }
          break;
        case AST_FUNCTION_ARCSEC:
          /* TRACE(("operate arcsec\n")); */
          if(1.0/stack[pos-1] > 1){
            stack[pos-1] = acos(1);
          }else if(1.0/stack[pos-1] < -1){
            stack[pos-1] = acos(-1);
          }else{
            stack[pos-1] = acos(1.0/stack[pos-1]);
          }
          break;
        case AST_FUNCTION_ARCCOT:
          /* TRACE(("operate arccot\n")); */
          if(eq->operator[i-1] == AST_MINUS
              && DOUBLE_EQ(*eq->number[i-2], 0)
              && DOUBLE_EQ(*eq->number[i-3], 0)){
            stack[pos-1] = atan(-1.0/stack[pos-1]);
          }else{
            stack[pos-1] = atan(1.0/stack[pos-1]);
          }
          break;
        case AST_FUNCTION_SINH:
          /* TRACE(("operate sinh\n")); */
          stack[pos-1] = sinh(stack[pos-1]);
          break;
        case AST_FUNCTION_COSH:
          /* TRACE(("operate cosh\n")); */
          stack[pos-1] = cosh(stack[pos-1]);
          break;
        case AST_FUNCTION_TANH:
          /* TRACE(("operate tanh\n")); */
          stack[pos-1] = tanh(stack[pos-1]);
          break;
        case AST_FUNCTION_CSCH:
          /* TRACE(("operate csch\n")); */
          stack[pos-1] = sinh(1.0/stack[pos-1]);
          break;
        case AST_FUNCTION_SECH:
          /* TRACE(("operate sech\n")); */
          stack[pos-1] = cosh(1.0/stack[pos-1]);
          break;
        case AST_FUNCTION_COTH:
          /* TRACE(("operate coth\n")); */
          stack[pos-1] = tanh(1.0/stack[pos-1]);
          break;
        case AST_FUNCTION_ARCSINH:
          /* TRACE(("operate arcsinh\n")); */
          stack[pos-1] = asinh(stack[pos-1]);
          break;
        case AST_FUNCTION_ARCCOSH:
          /* TRACE(("operate arccosh\n")); */
          stack[pos-1] = acosh(stack[pos-1]);
          break;
        case AST_FUNCTION_ARCTANH:
          /* TRACE(("operate arctanh\n")); */
          if(stack[pos-1] >= 1){
            stack[pos-1] = DBL_MAX;
          }else if(stack[pos-1] <= -1){
            stack[pos-1] = -DBL_MAX;
          }else{
            stack[pos-1] = atanh(stack[pos-1]);
          }
          break;
        case AST_FUNCTION_ARCCSCH:
          /* TRACE(("operate arccsch\n")); */
          stack[pos-1] = asinh(1.0/stack[pos-1]);
          break;
        case AST_FUNCTION_ARCSECH:
          /* TRACE(("operate arcsech\n")); */
          if(DOUBLE_EQ(stack[pos-1], 0)){
            stack[pos-1] = DBL_MAX;
          }else if(stack[pos-1] > 1){
            stack[pos-1] = 0;
          }else{
            stack[pos-1] = acosh(1.0/stack[pos-1]);
          }
          break;
        case AST_FUNCTION_ARCCOTH:
          /* TRACE(("operate arccoth\n")); */
          if(1.0/stack[pos-1] >= 1){
            stack[pos-1] = DBL_MAX;
          }else if(1.0/stack[pos-1] <= -1){
            stack[pos-1] = -DBL_MAX;
          }else{
            stack[pos-1] = atanh(1.0/stack[pos-1]);
          }
          break;
        case AST_FUNCTION_EXP:
          /* TRACE(("operate exp\n")); */
          stack[pos-1] = exp(stack[pos-1]);
          break;
        case AST_FUNCTION_LN:
          /* TRACE(("operate ln\n")); */
          stack[pos-1] = log(stack[pos-1]);
          break;
        case AST_FUNCTION_LOG:
          /* TRACE(("operate log\n")); */
          stack[pos-2] = log(stack[pos-1])/log(stack[pos-2]);
          pos--;
          break;
        case AST_FUNCTION_ROOT:
          /* TRACE(("operate root\n")); */
          stack[pos-2] = pow(stack[pos-1], 1/stack[pos-2]);
          pos--;
          break;
        case AST_FUNCTION_CEILING:
          /* TRACE(("operate ceiling\n")); */
          stack[pos-1] = ceil(stack[pos-1]);
          break;
        case AST_FUNCTION_FLOOR:
          /* TRACE(("operate floor\n")); */
          stack[pos-1] = floor(stack[pos-1]);
          break;
        case AST_FUNCTION_DELAY:
          /* TRACE(("operate delay\n")); */
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
          /* TRACE(("operate eq\n")); */
          /* TRACE(("EQ stack[pos-2] = %lf : stack[pos-1] = %lf\n", stack[pos-2], stack[pos-1])); */
          if(DOUBLE_EQ(stack[pos-2], stack[pos-1])){
            stack[pos-2] = 1;
          }else{
            stack[pos-2] = 0;
          }
          pos--;
          break;
        case AST_RELATIONAL_NEQ:
          /* TRACE(("operate neq\n")); */
          /* TRACE(("NEQ stack[pos-2] = %lf : stack[pos-1] = %lf\n", stack[pos-2], stack[pos-1])); */
          if(!DOUBLE_EQ(stack[pos-2], stack[pos-1])){
            stack[pos-2] = 1;
          }else{
            stack[pos-2] = 0;
          }
          pos--;
          break;
        case AST_RELATIONAL_LT:
          /* TRACE(("operate lt\n")); */
          if(stack[pos-2] < stack[pos-1]){
            stack[pos-2] = 1;
          }else{
            stack[pos-2] = 0;
          }
          pos--;
          break;
        case AST_RELATIONAL_GT:
          /* TRACE(("operate gt\n")); */
          if(stack[pos-2] > stack[pos-1]){
            stack[pos-2] = 1;
          }else{
            stack[pos-2] = 0;
          }
          pos--;
          break;
        case AST_RELATIONAL_LEQ:
          /* TRACE(("operate leq\n")); */
          if(stack[pos-2] <= stack[pos-1]){
            stack[pos-2] = 1;
          }else{
            stack[pos-2] = 0;
          }
          pos--;
          break;
        case AST_RELATIONAL_GEQ:
          /* TRACE(("operate geq\n")); */
          if(stack[pos-2] >= stack[pos-1]){
            stack[pos-2] = 1;
          }else{
            stack[pos-2] = 0;
          }
          pos--;
          break;
        case AST_LOGICAL_AND:
          /* TRACE(("operate and\n")); */
          if(stack[pos-2] >= 0.5 && stack[pos-1] >= 0.5){
            stack[pos-2] = 1;
          }else{
            stack[pos-2] = 0;
          }
          pos--;
          break;
        case AST_LOGICAL_NOT:
          /* TRACE(("operate not\n")); */
          if(stack[pos-1] >= 0.5){
            stack[pos-1] = 0;
          }else{
            stack[pos-1] = 1;
          }
          break;
        case AST_LOGICAL_OR:
          /* TRACE(("operate or\n")); */
          if(stack[pos-2] >= 0.5 || stack[pos-1] >= 0.5){
            stack[pos-2] = 1;
          }else{
            stack[pos-2] = 0;
          }
          pos--;
          break;
        case AST_LOGICAL_XOR:
          /* TRACE(("operate xor\n")); */
          if((stack[pos-2] >= 0.5 && stack[pos-1] < 0.5)
              || (stack[pos-2] < 0.5 && stack[pos-1] >= 0.5)){
            stack[pos-2] = 1;
          }else{
            stack[pos-2] = 0;
          }
          pos--;
          break;
        case AST_CONSTANT_TRUE:
          /* TRACE(("stack true\n")); */
          stack[pos] = 1;
          pos++;
          break;
        case AST_CONSTANT_FALSE:
          /* TRACE(("stack false\n")); */
          stack[pos] = 0;
          pos++;
          break;
      }
    }
  }

  return stack[0];
}
