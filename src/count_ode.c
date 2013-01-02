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
#include <time.h>
#include "libsbmlsim/libsbmlsim.h"

int count_ode(mySpecies* sp[], unsigned int num_of_species, int* ode_check, Species_t* s){
	unsigned int i;
	for(i=0; i<num_of_species; i++){
		if(strcmp(Species_getId(s), Species_getId(sp[i]->origin)) == 0){
			if(*(ode_check + i) == 0){
				*(ode_check + i) = 1;
				return 1;
			}
		}
	}
	return 0;
}
