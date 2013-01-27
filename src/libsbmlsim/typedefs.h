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
#ifndef LibSBMLSim_Typedefs_h
#define LibSBMLSim_Typedefs_h

typedef struct _equation equation;
typedef struct _mySpecies mySpecies;
typedef struct _myCompartment myCompartment;
typedef struct _myParameter myParameter;
typedef struct _myRule myRule;
typedef struct _mySpeciesReference mySpeciesReference;
typedef struct _myReaction myReaction;
typedef struct _myEvent myEvent;
typedef struct _myEventAssignment myEventAssignment;
typedef struct _myDelay myDelay;
typedef struct _myInitialAssignment myInitialAssignment;
typedef struct _allocated_memory allocated_memory;
typedef struct _copied_AST copied_AST;

/* no header files yet */
typedef struct _timeVariantAssignments timeVariantAssignments;
typedef struct _myASTNode myASTNode;
typedef struct _myAlgebraicEquations myAlgebraicEquations;
typedef struct _myAlgTargetSp myAlgTargetSp;
typedef struct _myAlgTargetParam myAlgTargetParam;
typedef struct _myAlgTargetComp myAlgTargetComp;

#endif  /* LibSBMLSim_Typedefs_h */
