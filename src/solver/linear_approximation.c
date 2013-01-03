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

double approximate_delay_linearly(double* stack, int pos, double** delay_preserver, double* time, int rk_order, myResult* res, int cycle, int print_interval, int* err_zero_flag) {
	int i;
	double grad = 0.0;
	double result_value = 0.0;
	double delayed_time = *(time) - stack[pos - 1];
	double* values_time = res->values_time_fordelay;
	if (*(err_zero_flag) == 0) {
		/* calclulate gradient -> linear approximation */
		for (i=1; i<cycle; i++) {
			if ( *(values_time + i - 1) <= delayed_time && delayed_time < *(values_time + i )) {
				grad =  ( *(*(delay_preserver + i) + rk_order) - *( *(delay_preserver + (i - 1)) + rk_order)) / (*(values_time + i) - *(values_time + (i-1)));
				result_value = *(*(delay_preserver + (i - 1)) + rk_order )  + grad * (delayed_time - *(values_time + i-1));
				break;
			}
		}
	}else {
		/* calclulate gradient -> linear approximation */
		for (i=1; i<cycle; i++) {
			if ( *(values_time + i - 1) <= delayed_time && delayed_time < *(values_time + i )) {
				grad =  ( *(*(delay_preserver + i * print_interval) + rk_order) - *( *(delay_preserver + (i - 1) * print_interval) + rk_order)) / (*(values_time + i) - *(values_time + (i-1)));
				result_value = *(*(delay_preserver + (i - 1) * print_interval) + rk_order )  + grad * (delayed_time - *(values_time + i-1));
				break;
			}
		}
	}
	return result_value;
}


double approximate_printresult_linearly(double value, double temp_value, double value_time, double tempvalue_time, double fixed_time) {
	double grad = 0.0;
	double result = 0.0;
	grad = (temp_value - value) / (value_time - tempvalue_time);
	result = value + grad * (fixed_time - tempvalue_time);
	return result;

}
