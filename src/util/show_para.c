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

void show_para(Model_t *m){
		  unsigned int i, j;
  int column = 0;
  printf("(Parameter List)\n");
  for(i=0; i<Model_getNumParameters(m); i++){
	  printf("[%d]global : ID=%s Name=%s value=%.16g\n", column, Parameter_getId((Parameter_t*)ListOf_get(Model_getListOfParameters(m), i)), Parameter_getName((Parameter_t*)ListOf_get(Model_getListOfParameters(m), i)), Parameter_getValue((Parameter_t*)ListOf_get(Model_getListOfParameters(m), i)));
	  column++;
  }
  for(i = 0; i < Model_getNumReactions(m); i++) {
	  for(j = 0; j < ListOf_size(KineticLaw_getListOfParameters(Reaction_getKineticLaw((Reaction_t*)ListOf_get(Model_getListOfReactions(m), i)))); j++){
		  printf("[%d]local  : ID=%s Name=%s value=%.16g\n", column, Parameter_getId(KineticLaw_getParameter(Reaction_getKineticLaw((Reaction_t*)ListOf_get(Model_getListOfReactions(m), i)), j)), Parameter_getName(KineticLaw_getParameter(Reaction_getKineticLaw((Reaction_t*)ListOf_get(Model_getListOfReactions(m), i)), j)), Parameter_getValue(KineticLaw_getParameter(Reaction_getKineticLaw((Reaction_t*)ListOf_get(Model_getListOfReactions(m), i)), j)));
		  column++;
	  }
  }
}
