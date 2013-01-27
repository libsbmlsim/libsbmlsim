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
#ifndef LibSBMLSim_Common_h
#define LibSBMLSim_Common_h

/* debug print */
#ifdef DEBUG_PRINT
#define DEBUG_PRINT_FLAG 1
#else
#define DEBUG_PRINT_FLAG 0
#endif

#define TRACE(x) do { if (DEBUG_PRINT_FLAG) dbg_printf x; } while (0)

/* progress print */
#ifdef PROGRESS_PRINT
#define PROGRESS_PRINT_FLAG 1
#else
#define PROGRESS_PRINT_FLAG 0
#endif

#define PRG_TRACE(x) do { if (PROGRESS_PRINT_FLAG) prg_printf x; } while (0)

/* defined variables */
#define MAX_MATH_LENGTH 256
#define MAX_DELAY_REACTION_NUM 256
#define MAX_ARG_NUM 256
#define MAX_ALLOCATED_MEMORY 1024
#define MAX_COPIED_AST 1024
#define MAX_ALGEBRAIC_VARIABLES 1024
#define MAX_ALGEBRAIC_CONSTANTS 1024
#define MAX_IDENTICAL_EVENTS 256
#define MAX_EVENTASSIGNMENTS 256
#define MAX_TIME_VARIANT_ASSIGNMENT 1024
#define MAX_INCLUDING_SPECIES 256
#define MAX_MTHD_NAME_LENGTH 256

#define EPSIRON 1.0e-8
#define ABSOLUTE_ERROR_TOLERANCE 1.0e-9
#define RELATIVE_ERROR_TOLERANCE 1.0e-6
#define FINE_ABSOLUTE_ERROR_TOLERANCE 1.0e-22
#define FINE_RELATIVE_ERROR_TOLERANCE 1.0e-11
#define DEFAULT_FACMAX 2.0
#define DOUBLE_EQ(x, v) (((v - EPSIRON) < x) && (x < (v + EPSIRON)))

#endif  /* LibSBMLSim_Common_h */
