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
#include "../libsbmlsim/libsbmlsim.h"

static int list_has_element(char *list[], int list_length, char *element){
  int i;
  int flag = 0;
  for(i=0; i<list_length; i++){
    if(strcmp(element, list[i]) == 0){
      flag = 1;
    }
  }
  return flag;
}

int assign_ok(ASTNode_t *assignment_math, char *target_list[], int num_of_targets, char* assigned_target_list[], unsigned int num_of_assigned_targets, int flag){
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

void calc_initial_assignment(myInitialAssignment *initAssign[], unsigned int num_of_initialAssignments, double dt, int cycle, double *reverse_time){
  unsigned int i;
  char **target_list;
  char **assigned_target_list;
  ASTNode_t **assignment_math_list;
  unsigned int num_of_assigned_targets = 0;

  target_list = (char **)malloc(sizeof(char *) * num_of_initialAssignments);
  assigned_target_list = (char **)malloc(sizeof(char *) * num_of_initialAssignments);
  assignment_math_list = (ASTNode_t **)malloc(sizeof(ASTNode_t *) * num_of_initialAssignments);

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
  free(target_list);
  free(assigned_target_list);
  free(assignment_math_list);
}

void calc_initial_assignmentf(myInitialAssignment *initAssign[], unsigned int num_of_initialAssignments, double dt, int cycle, double *reverse_time, double* time, myResult* result, int print_interval, int* err_zero_flag){
  unsigned int i;
  char **target_list;
  char **assigned_target_list;
  ASTNode_t **assignment_math_list;
  unsigned int num_of_assigned_targets = 0;
  target_list = (char **)malloc(sizeof(char *) * num_of_initialAssignments);
  assigned_target_list = (char **)malloc(sizeof(char *) * num_of_initialAssignments);
  assignment_math_list = (ASTNode_t **)malloc(sizeof(ASTNode_t *) * num_of_initialAssignments);

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
			  initAssign[i]->target_species->temp_value = calcf(initAssign[i]->eq, dt, cycle, reverse_time, 0, time, time, result, print_interval, err_zero_flag);
            initAssign[i]->target_species->value = initAssign[i]->target_species->temp_value;
          }else if(initAssign[i]->target_parameter != NULL){
			  initAssign[i]->target_parameter->temp_value = calcf(initAssign[i]->eq, dt, cycle, reverse_time, 0, time, time, result, print_interval, err_zero_flag);
            initAssign[i]->target_parameter->value = initAssign[i]->target_parameter->temp_value;
          }else if(initAssign[i]->target_compartment != NULL){
			  initAssign[i]->target_compartment->temp_value = calcf(initAssign[i]->eq, dt, cycle, reverse_time, 0, time, time, result, print_interval, err_zero_flag);
            initAssign[i]->target_compartment->value = initAssign[i]->target_compartment->temp_value;
          }else if(initAssign[i]->target_species_reference != NULL){
			  initAssign[i]->target_species_reference->temp_value = calcf(initAssign[i]->eq, dt, cycle, reverse_time, 0, time, time, result, print_interval, err_zero_flag);
            initAssign[i]->target_species_reference->value = initAssign[i]->target_species_reference->temp_value;
          }
          assigned_target_list[num_of_assigned_targets++] = target_list[i];
          TRACE(("target : %s is assigned to %lf\n", target_list[i], calcf(initAssign[i]->eq, dt, cycle, reverse_time, 0, time, time, result, print_interval, err_zero_flag)));
        }
      }
    }
  }
  free(target_list);
  free(assigned_target_list);
  free(assignment_math_list);
}

