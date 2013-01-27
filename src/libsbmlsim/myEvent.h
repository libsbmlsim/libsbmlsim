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
#ifndef LibSBMLSim_MyEvent_h
#define LibSBMLSim_MyEvent_h

#include "typedefs.h"
#include "boolean.h"
#include "equation.h"
#include "myEventAssignment.h"
#include "myDelay.h"
#include <sbml/SBMLTypes.h>

struct _myEvent {
  Event_t *origin;
  equation *eq; /* condition equation */
  myEventAssignment** assignments;
  boolean is_able_to_fire;
  myDelay *event_delay;
  double *firing_times;
  unsigned int num_of_delayed_events_que;
  int next_firing_index;
  boolean is_persistent;
  equation *priority_eq;
};

myEvent *myEvent_create();
void myEvent_initWithModel(myEvent *event, Model_t *mode, int index);
void myEvent_free(myEvent *event);

Event_t *myEvent_getOrigin(myEvent *event);

#endif /* LibSBMLSim_MyEvent_h */
