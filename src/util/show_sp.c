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

void show_sp(Model_t *m){
  unsigned int i;
  int column = 0;
  printf("(Species List)\n");
  for(i=0; i<Model_getNumSpecies(m); i++){
	  printf("column %d : ID=%s Name=%s Initial Amount=%.16g\n", column, Species_getId((Species_t*)ListOf_get(Model_getListOfSpecies(m), i)), Species_getName((Species_t*)ListOf_get(Model_getListOfSpecies(m), i)), Species_getInitialAmount((Species_t*)ListOf_get(Model_getListOfSpecies(m), i)));
	  column++;
  }
}

