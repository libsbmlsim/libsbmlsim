#include "libsbmlsim/libsbmlsim.h"

int list_has_element(char *list[], int list_length, char *element){
  int i;
  int flag = 0;
  for(i=0; i<list_length; i++){
    if(strcmp(element, list[i]) == 0){
      flag = 1;
    }
  }
  return flag;
}

int assign_ok(ASTNode_t *assignment_math, char *target_list[], int num_of_targets, char* assigned_target_list[], int num_of_assigned_targets, int flag){
  ASTNode_t *left, *right;
  char *name;
  if((left=ASTNode_getLeftChild(assignment_math)) != NULL){
    flag = assign_ok(left, target_list, num_of_targets, assigned_target_list, num_of_assigned_targets, flag);
  }
  if((right=ASTNode_getRightChild(assignment_math)) != NULL){
    flag = assign_ok(right, target_list, num_of_targets, assigned_target_list, num_of_assigned_targets, flag);
  }
  if(ASTNode_getType(assignment_math) == AST_NAME){
    name = (char*)ASTNode_getName(assignment_math);
    if(list_has_element(target_list, num_of_targets, name)){
      if(!list_has_element(assigned_target_list, num_of_assigned_targets, name)){
        flag = 0;
      }
    }
  }
  return flag;
}

void calc_initial_assignment(myInitialAssignment *initAssign[], int num_of_initialAssignments, double dt, int cycle, double *reverse_time){
  int i;
  char *target_list[num_of_initialAssignments];
  char *assigned_target_list[num_of_initialAssignments];
  int num_of_assigned_targets = 0;
  ASTNode_t *assignment_math_list[num_of_initialAssignments];

  for(i=0; i<num_of_initialAssignments; i++){
    if(initAssign[i]->target_species != NULL){
      target_list[i] = (char*)Species_getId(initAssign[i]->target_species->origin);
    }
    if(initAssign[i]->target_parameter != NULL){
      target_list[i] = (char*)Parameter_getId(initAssign[i]->target_parameter->origin);
    }
    if(initAssign[i]->target_compartment != NULL){
      target_list[i] = (char*)Compartment_getId(initAssign[i]->target_compartment->origin);
    }
    if(initAssign[i]->target_species_reference != NULL){
      target_list[i] = (char*)SpeciesReference_getId(initAssign[i]->target_species_reference->origin);
    }
  }

  for(i=0; i<num_of_initialAssignments; i++){
    assignment_math_list[i] = (ASTNode_t*)InitialAssignment_getMath(initAssign[i]->origin);
  }

  while(num_of_assigned_targets < num_of_initialAssignments){
    for(i=0; i<num_of_initialAssignments; i++){
      if(!list_has_element(assigned_target_list, num_of_assigned_targets, target_list[i])){
        if(assign_ok(assignment_math_list[i], target_list, num_of_initialAssignments, assigned_target_list, num_of_assigned_targets, 1)){
          if(initAssign[i]->target_species != NULL){
            initAssign[i]->target_species->temp_value = calc(initAssign[i]->eq, dt, cycle, reverse_time, 0);
            initAssign[i]->target_species->value = initAssign[i]->target_species->temp_value;
          }else if(initAssign[i]->target_parameter != NULL){
            initAssign[i]->target_parameter->temp_value = calc(initAssign[i]->eq, dt, cycle, reverse_time, 0);
            initAssign[i]->target_parameter->value = initAssign[i]->target_parameter->temp_value;
          }else if(initAssign[i]->target_compartment != NULL){
            initAssign[i]->target_compartment->temp_value = calc(initAssign[i]->eq, dt, cycle, reverse_time, 0);
            initAssign[i]->target_compartment->value = initAssign[i]->target_compartment->temp_value;
          }else if(initAssign[i]->target_species_reference != NULL){
            initAssign[i]->target_species_reference->temp_value = calc(initAssign[i]->eq, dt, cycle, reverse_time, 0);
            initAssign[i]->target_species_reference->value = initAssign[i]->target_species_reference->temp_value;
          }
          assigned_target_list[num_of_assigned_targets++] = target_list[i];
          TRACE(("target : %s is assigned to %lf\n", target_list[i], calc(initAssign[i]->eq, dt, cycle, reverse_time, 0)));
        }
      }
    }
  }
}
