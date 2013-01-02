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
#include "libsbmlsim/libsbmlsim.h"

void print_result_list(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[]){
  FILE *fp = NULL;
  unsigned int i;
  int column = 2;
  if ((fp = my_fopen(fp, "./simulation_results/result_list.dat", "w")) == NULL) {
    return;
  }
  TRACE(("Result : Species List\n"));
  fprintf(fp, "Result : Species List\n");
  for(i=0; i<Model_getNumSpecies(m); i++){
    /*     if(!(Species_getConstant(mySp[i]->origin) && Species_getBoundaryCondition(mySp[i]->origin))){ //XXX must remove this */
    fprintf(fp, "column %d : ID=%s Name=%s\n", column, Species_getId(mySp[i]->origin), Species_getName(mySp[i]->origin));
    column++;
    /*    } */
  }
  TRACE(("\n"));
  fprintf(fp, "\n");
  fprintf(fp, "Result : Parameter List\n");
  column = 2;
  for(i=0; i<Model_getNumParameters(m); i++){
    /*     if(!Parameter_getConstant(myParam[i]->origin)){ //XXX must remove this */
    TRACE(("column %d : ID=%s Name=%s\n", column, Parameter_getId(myParam[i]->origin), Parameter_getName(myParam[i]->origin)));
    fprintf(fp, "column %d : ID=%s Name=%s\n", column, Parameter_getId(myParam[i]->origin), Parameter_getName(myParam[i]->origin));
    column++;
    /*     } */
  }
  TRACE(("Result : Compartment List\n"));
  fprintf(fp, "Result : Compartment List\n");
  column = 2;
  for(i=0; i<Model_getNumCompartments(m); i++){
    /*     if(!Compartment_getConstant(myComp[i]->origin)){ //XXX must remove this */
    TRACE(("column %d : ID=%s Name=%s\n", column, Compartment_getId(myComp[i]->origin), Compartment_getName(myComp[i]->origin)));
    fprintf(fp, "column %d : ID=%s Name=%s\n", column, Compartment_getId(myComp[i]->origin), Compartment_getName(myComp[i]->origin));
    column++;
    /*     } */
  }
  fclose(fp);
}
