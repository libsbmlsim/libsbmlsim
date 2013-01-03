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

void forwarding_value(mySpecies *sp[], int sp_num, myParameter *param[], int param_num, myCompartment *comp[], int comp_num, mySpeciesReference *spr[], int spr_num){
  int i;
  for(i=0; i<sp_num; i++){
    sp[i]->value = sp[i]->temp_value;
  }
  for(i=0; i<param_num; i++){
    param[i]->value = param[i]->temp_value;
  }
  for(i=0; i<comp_num; i++){
    comp[i]->value = comp[i]->temp_value;
  }
  for(i=0; i<spr_num; i++){
    spr[i]->value = spr[i]->temp_value;
  }  
}
