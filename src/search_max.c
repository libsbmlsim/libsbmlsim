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

double search_max(myResult* result, int sta_var_column){
	int i;
	double *value_time_p = result->values_time;
	double *value_sp_p = result->values_sp;
	double max = DBL_MIN;
	for(i = 0; i < result->num_of_rows; i++) {
		if (i == 0) {
			max = *(value_sp_p + sta_var_column + result->num_of_columns_sp * i);
		}
		value_time_p++;
		/* calculate the maximum of state_variable */
		if (max < *(value_sp_p + sta_var_column + result->num_of_columns_sp * i)) {
			max = *(value_sp_p + sta_var_column + result->num_of_columns_sp * i);
		}
	}
	return max;
}

double search_local_max(myResult* result, int sta_var_column, double transition_time, double sim_time){
	int i;
	double *value_time_p = result->values_time;
	double *value_sp_p = result->values_sp;
	double local_max = DBL_MIN;
	for(i = 0; i < result->num_of_rows; i++) {
		if (transition_time < *(value_time_p) && *(value_time_p) < sim_time && transition_time > *(value_time_p - 1)) {
			local_max = *(value_sp_p + sta_var_column + result->num_of_columns_sp * i);
		}
		value_time_p++;
		/* calculate the local maximum of state_variable */
		if (local_max < *(value_sp_p + sta_var_column + result->num_of_columns_sp * i) && transition_time < *(value_time_p) && *(value_time_p) < sim_time) {
			local_max = *(value_sp_p + sta_var_column + result->num_of_columns_sp * i);
		}
	}
	return local_max;
}

double search_local_min(myResult* result, int sta_var_column, double transition_time, double sim_time){
	int i;
	double *value_time_p = result->values_time;
	double *value_sp_p = result->values_sp;
	double local_min = DBL_MAX;
	for(i = 0; i < result->num_of_rows; i++) {
		if (transition_time < *(value_time_p) && *(value_time_p) < sim_time && transition_time > *(value_time_p - 1)) {
			local_min = *(value_sp_p + sta_var_column + result->num_of_columns_sp * i);
		}
		value_time_p++;
		/* calculate the local minimum of state_variable */
		if (local_min > *(value_sp_p + sta_var_column + result->num_of_columns_sp * i) && transition_time < *(value_time_p) && *(value_time_p) < sim_time) {
			local_min = *(value_sp_p + sta_var_column + result->num_of_columns_sp * i);
		}
	}
	return local_min;
}
