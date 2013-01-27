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
#include "libsbmlsim/myRule.h"
#include <stdlib.h>
#include <sbml/SBMLTypes.h>

myRule *myRule_create() {
  myRule *rule = (myRule *)malloc(sizeof(myRule));
  rule->origin = NULL;
  rule->eq = NULL;
  rule->target_species = NULL;
  rule->target_parameter = NULL;
  rule->target_compartment = NULL;
  rule->target_species_reference = NULL;
  rule->is_rate = false;
  rule->is_assignment = false;
  rule->is_algebraic = false;
  return rule;
}

void myRule_initWithModel(myRule *rule, Model_t *model, int index) {
  Rule_t *origin = (Rule_t *)ListOf_get(Model_getListOfRules(model), index);

  rule->origin = origin;
  rule->eq = equation_create();
  rule->is_rate = Rule_isRate(origin);
  rule->is_assignment = Rule_isAssignment(origin);
  rule->is_algebraic = Rule_isAlgebraic(origin);
}

void myRule_free(myRule *rule) {
  if (rule == NULL) {
    return;
  }

  /*
  if (!rule->is_algebraic) {
    equation_free(rule->eq);
  }
  */
  if (rule->eq != NULL) {
    equation_free(rule->eq);
  }
  free(rule);
}

void myRule_initTarget(
    myRule *rule,
    mySpecies **species, unsigned int num_of_species,
    myParameter **parameters, unsigned int num_of_parameters,
    myCompartment **compartments, unsigned int num_of_compartments,
    myReaction **reactions, unsigned int num_of_reactions) {
  unsigned int i, j;
  Rule_t *origin;
  const char *origin_var;
  mySpecies *sp;
  myParameter *param;
  myCompartment *comp;
  myReaction *re;
  mySpeciesReference *product, *reactant;

  origin = myRule_getOrigin(rule);
  origin_var = Rule_getVariable(origin);

  for (i = 0; i < num_of_species; i++) {
    sp = species[i];
    if (strcmp(origin_var, Species_getId(sp->origin)) == 0) {
      myRule_setTargetSpecies(rule, sp);
      mySpecies_setDependingRule(sp, rule);
      return;
    }
  }

  for (i = 0; i < num_of_parameters; i++) {
    param = parameters[i];
    if (strcmp(origin_var, Parameter_getId(param->origin)) == 0) {
      myRule_setTargetParameter(rule, param);
      myParameter_setDependingRule(param, rule);
      return;
    }
  }

  for (i = 0; i < num_of_compartments; i++) {
    comp = compartments[i];
    if (strcmp(origin_var, Compartment_getId(comp->origin)) == 0) {
      myRule_setTargetCompartment(rule, comp);
      myCompartment_setDependingRule(comp, rule);
      return;
    }
  }

  for (i=0; i < num_of_reactions; i++) {
    re = reactions[i];

    for (j = 0; j < re->num_of_products; j++) {
      product = re->products[j];
      if (SpeciesReference_isSetId(product->origin)
          && strcmp(origin_var, SpeciesReference_getId(product->origin)) == 0) {
        myRule_setTargetSpeciesReference(rule, product);
        mySpeciesReference_setDependingRule(product, rule);
        return;
      }
    }

    for (j = 0; j < re->num_of_reactants; j++) {
      reactant = re->reactants[j];
      if (SpeciesReference_isSetId(reactant->origin)
          && strcmp(origin_var, SpeciesReference_getId(reactant->origin)) == 0) {
        myRule_setTargetSpeciesReference(rule, reactant);
        mySpeciesReference_setDependingRule(reactant, rule);
        return;
      }
    }
  }
}

Rule_t *myRule_getOrigin(myRule *rule) {
  return rule->origin;
}

equation *myRule_getEquation(myRule *rule) {
  return rule->eq;
}

void myRule_setTargetSpecies(myRule *rule, mySpecies *species) {
  rule->target_species = species;
}

void myRule_setTargetParameter(myRule *rule, myParameter *parameter) {
  rule->target_parameter = parameter;
}

void myRule_setTargetCompartment(myRule *rule, myCompartment *compartment) {
  rule->target_compartment = compartment;
}

void myRule_setTargetSpeciesReference(myRule *rule, mySpeciesReference *ref) {
  rule->target_species_reference = ref;
}

