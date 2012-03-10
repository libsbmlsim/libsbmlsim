#include "sbml/common/common.h"
#include "sbml/SBMLReader.h"
#include "sbml/SBMLTypes.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include <float.h>
#include "header.h"

void recursive_calc_event(myEvent *event[], int num_of_events, myEvent *event_buf[], int *num_of_remained_events, double *assignment_values_from_trigger_time[], double dt, double time, int cycle, double *reverse_time){
  int i, j, k;
  int is_condition_satisfied;
  int flag;
  myEvent *temp_event;
  double priority_buf[MAX_IDENTICAL_EVENTS];
  double temp_priority;
  int num_of_same_priority_events;
  int selected_order;
  double *temp_assignment_values_from_trigger_time;

  //event_buf check(remove persistent = false && trigger = false)
  for(i=0; i<(*num_of_remained_events); i++){
    if(!event_buf[i]->is_persistent && calc(event_buf[i]->eq, dt, cycle, reverse_time, 0) < 0.5){
      //dbg_printf("%s is deleted\n", Event_getId(event_buf[i]->origin));
      for(j=i; j<(*num_of_remained_events)-1; j++){
	for(k=0; k<Event_getNumEventAssignments(event_buf[j+1]->origin); k++){
	  assignment_values_from_trigger_time[j][k] = assignment_values_from_trigger_time[j+1][k];
	}
	event_buf[j] = event_buf[j+1];
      }
      (*num_of_remained_events)--;
    }
  }
  //event cycle
  for(i=0; i<num_of_events; i++){
    flag = 0;
    //condition determinatioin
    if(calc(event[i]->eq, dt, cycle, reverse_time, 0) >= 0.5){
      is_condition_satisfied = 1;
    }else{
      is_condition_satisfied = 0;
      event[i]->is_able_to_fire = 1;
    }
    //firing flag determination
    if(is_condition_satisfied && event[i]->is_able_to_fire){
      if(event[i]->event_delay == NULL){
	flag = 1;
	event[i]->is_able_to_fire = 0;
      }else{
	event[i]->firing_times[event[i]->num_of_delayed_events_que++] = time + calc(event[i]->event_delay->eq, dt, cycle, reverse_time, 0);
	event[i]->is_able_to_fire = 0;
      }
    }
    //delay event firing determination
    if(event[i]->event_delay != NULL && time >= event[i]->firing_times[event[i]->next_firing_index]){
      if(event[i]->is_persistent || calc(event[i]->eq, dt, cycle, reverse_time, 0) >= 0.5){
	flag = 1;
      }
      event[i]->next_firing_index++;
    }
    if(flag){
      //buffering
      //dbg_printf("%s is buffered in event_buf[%d]\n", Event_getId(event[i]->origin), *num_of_remained_events);
      event_buf[(*num_of_remained_events)] = event[i];
      for(j=0; j<Event_getNumEventAssignments(event[i]->origin); j++){
	assignment_values_from_trigger_time[(*num_of_remained_events)][j] = calc(event[i]->assignments[j]->eq, dt, cycle, reverse_time, 0);
	//dbg_printf("%lf(%p) is buffered in assignment_values_from_trigger_time[%d][%d]\n", assignment_values_from_trigger_time[(*num_of_remained_events)][j], &(assignment_values_from_trigger_time[(*num_of_remained_events)][j]), *num_of_remained_events, j);
      }
      (*num_of_remained_events)++;
    }
  }
  //calculate priority
  for(i=0; i<(*num_of_remained_events); i++){
    if(event_buf[i]->priority_eq != NULL){
      priority_buf[i] = calc(event_buf[i]->priority_eq, dt, cycle, reverse_time, 0);
    }else{
      priority_buf[i] = -DBL_MAX;
    }
  }
  //sort
  for(i=0; i<(*num_of_remained_events); i++){
    for(j=(*num_of_remained_events)-1; j>i; j--){
      if(priority_buf[j] >= priority_buf[j-1]){
	//swap priority_buf
	temp_priority = priority_buf[j-1];
	priority_buf[j-1] = priority_buf[j];
	priority_buf[j] = temp_priority;
	//swap event_buf
	temp_event = event_buf[j-1];
	event_buf[j-1] = event_buf[j];
	event_buf[j] = temp_event;
	//swap assignment_values_from_trigger_time
	temp_assignment_values_from_trigger_time = assignment_values_from_trigger_time[j-1];
	assignment_values_from_trigger_time[j-1] = assignment_values_from_trigger_time[j];
	assignment_values_from_trigger_time[j] = temp_assignment_values_from_trigger_time;
      }
    }
  }
  //count same priority events
  num_of_same_priority_events = 1;
  for(i=1; i<(*num_of_remained_events); i++){
    if(priority_buf[i-1] == priority_buf[i]){
      num_of_same_priority_events++;
    }else{
      break;
    }
  }
  //rondom order determination for same priority events
  selected_order = rand()%num_of_same_priority_events;
  //swap event_buf
  temp_event = event_buf[0];
  event_buf[0] = event_buf[selected_order];
  event_buf[selected_order] = temp_event;
  //swap assignment
  temp_assignment_values_from_trigger_time = assignment_values_from_trigger_time[0];
  assignment_values_from_trigger_time[0] = assignment_values_from_trigger_time[selected_order];
  assignment_values_from_trigger_time[selected_order] = temp_assignment_values_from_trigger_time;
}

void calc_event(myEvent *event[], int num_of_events, double dt, double time, int cycle, double *reverse_time){
  int i, j;
  myEvent *event_buf[MAX_IDENTICAL_EVENTS];
  static double sub_assignment_values_from_trigger_time[MAX_IDENTICAL_EVENTS][MAX_EVENTASSIGNMENTS];
  static double *assignment_values_from_trigger_time[MAX_IDENTICAL_EVENTS];
  int num_of_remained_events = 0;
  myEventAssignment* assignment;

  if(cycle == 0){
    for(i=0; i<MAX_IDENTICAL_EVENTS; i++){
      assignment_values_from_trigger_time[i] = sub_assignment_values_from_trigger_time[i];
    }
  }

  //recursive processing
  recursive_calc_event(event, num_of_events, event_buf, &num_of_remained_events, assignment_values_from_trigger_time, dt, time, cycle, reverse_time);
 
  //proccess assignment start
  while(num_of_remained_events != 0){
    if(event_buf[0]->is_persistent || calc(event_buf[0]->eq, dt, cycle, reverse_time, 0) >= 0.5){
      //dbg_printf("%s's assignment is processed\n", Event_getId(event_buf[0]->origin));
      //forwarding value
      if(Event_getUseValuesFromTriggerTime(event_buf[0]->origin)){
	for(i=0; i<Event_getNumEventAssignments(event_buf[0]->origin); i++){
	  assignment = event_buf[0]->assignments[i];
	  //dbg_printf("value used in assignment is %lf(%p)\n", assignment_values_from_trigger_time[0][i], &assignment_values_from_trigger_time[0][i]);
	  if(assignment->target_species != NULL){
	    assignment->target_species->value = assignment_values_from_trigger_time[0][i];
	  }else if(assignment->target_parameter != NULL){
	    assignment->target_parameter->value = assignment_values_from_trigger_time[0][i];
	  }else if(assignment->target_compartment != NULL){
	    assignment->target_compartment->value = assignment_values_from_trigger_time[0][i];	       
	  }else if(assignment->target_species_reference != NULL){
	    assignment->target_species_reference->value = assignment_values_from_trigger_time[0][i];
	  }
	}
      }else{
	for(i=0; i<Event_getNumEventAssignments(event_buf[0]->origin); i++){
	  assignment = event_buf[0]->assignments[i];
	  if(assignment->target_species != NULL){
	    assignment->target_species->value = calc(assignment->eq, dt, cycle, reverse_time, 0);
	  }else if(assignment->target_parameter != NULL){
	    assignment->target_parameter->value = calc(assignment->eq, dt, cycle, reverse_time, 0);
	  }else if(assignment->target_compartment != NULL){
	    assignment->target_compartment->value = calc(assignment->eq, dt, cycle, reverse_time, 0);	       
	  }else if(assignment->target_species_reference != NULL){
	    assignment->target_species_reference->value = calc(assignment->eq, dt, cycle, reverse_time, 0);
	  }
	}
      }
      //forwarding temp value
      for(i=0; i<Event_getNumEventAssignments(event_buf[0]->origin); i++){
	assignment = event_buf[0]->assignments[i];
	if(assignment->target_species != NULL){
	  assignment->target_species->temp_value = assignment->target_species->value;
	}else if(assignment->target_parameter != NULL){
	  assignment->target_parameter->temp_value = assignment->target_parameter->value;
	}else if(assignment->target_compartment != NULL){
    //new code
    for(j=0; j<assignment->target_compartment->num_of_including_species; j++){
      if(assignment->target_compartment->including_species[j]->is_concentration){
        assignment->target_compartment->including_species[j]->value = assignment->target_compartment->including_species[j]->value*assignment->target_compartment->temp_value/assignment->target_compartment->value;
        assignment->target_compartment->including_species[j]->temp_value = assignment->target_compartment->including_species[j]->value;
      }
    }
    //
    assignment->target_compartment->temp_value = assignment->target_compartment->value;
	}else if(assignment->target_species_reference != NULL){
	  assignment->target_species_reference->temp_value = assignment->target_species_reference->value;
	}
      } 
    }
    for(i=1; i<num_of_remained_events; i++){
      for(j=0; j<Event_getNumEventAssignments(event_buf[i]->origin); j++){
	assignment_values_from_trigger_time[i-1][j] = assignment_values_from_trigger_time[i][j];
      }
      event_buf[i-1] = event_buf[i];
    }
    num_of_remained_events--;
/*     dbg_printf("buffer state is\n"); */
/*     for(i=0; i<num_of_remained_events; i++){ */
/*       dbg_printf("event_buf[%d] = %s\n", i, Event_getId(event_buf[i]->origin)); */
/*       for(j=0; j<Event_getNumEventAssignments(event_buf[i]->origin); j++){ */
/* 	dbg_printf("assignment_values_from_trigger_time[%d][%d] = %lf\n", i, j, assignment_values_from_trigger_time[i][j]); */
/*       } */
/*     } */
    //recursive processing
    recursive_calc_event(event, num_of_events, event_buf, &num_of_remained_events, assignment_values_from_trigger_time, dt, time, cycle, reverse_time);
  }//proccess assignment finish
    
}
