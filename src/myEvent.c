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
#include "libsbmlsim/myEvent.h"
#include <stdlib.h>
#include <sbml/SBMLTypes.h>

myEvent *myEvent_create() {
  myEvent *event = (myEvent *)malloc(sizeof(myEvent));
  event->origin = NULL;
  event->eq = NULL;
  event->assignments = NULL;
  event->is_able_to_fire = false;
  event->event_delay = NULL;
  event->firing_times = NULL;
  event->num_of_delayed_events_que = 0;
  event->next_firing_index = 0;
  event->is_persistent = false;
  event->priority_eq = NULL;
  return event;
}

void myEvent_initWithModel(myEvent *event, Model_t *model, int index) {
  Event_t *origin = (Event_t *)ListOf_get(Model_getListOfEvents(model), index);

  event->origin = origin;
  event->eq = equation_create();
  if (Trigger_isSetPersistent(Event_getTrigger(origin))) {
    event->is_persistent = Trigger_getPersistent(Event_getTrigger(origin));
  } else {
    event->is_persistent = true;
  }
  if (Trigger_isSetInitialValue(Event_getTrigger(origin))) {
    event->is_able_to_fire = !Trigger_getInitialValue(Event_getTrigger(origin));
  } else {
    event->is_able_to_fire = true;
  }
}

void myEvent_free(myEvent *event) {
  unsigned int i;

  if (event == NULL) {
    return;
  }

  if (event->assignments != NULL) {
    for (i = 0; i < Event_getNumEventAssignments(event->origin); i++) {
      myEventAssignment_free(event->assignments[i]);
    }
    free(event->assignments);
  }

  if (event->event_delay != NULL) {
    equation_free(event->event_delay->eq);
    free(event->event_delay);
  }

  if (event->firing_times != NULL) {
    free(event->firing_times);
  }

  if (event->priority_eq != NULL) {
    free(event->priority_eq);
  }

  if (event->eq != NULL) {
    free(event->eq);
  }

  free(event);
}

Event_t *myEvent_getOrigin(myEvent *event) {
  return event->origin;
}

