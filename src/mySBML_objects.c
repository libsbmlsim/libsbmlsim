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

/* private functions */
static int include_time(ASTNode_t *node, int flag);
static void locate_species_in_compartment(
    mySpecies **species, int num_of_species,
    myCompartment **compartments, int num_of_compartments);
static void create_species(mySpecies *species[], Model_t *model);
static void create_parameters(myParameter *parameters[], Model_t *model);
static void create_compartments(myCompartment *compartments[], Model_t *model);
static void create_reactions(myReaction *reactions[], mySpecies *species[], Model_t *model);
/*********************/


/* create my SBML obejects for efficient simulations */
void create_mySBML_objects(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], myReaction *myRe[], myRule *myRu[], myEvent *myEv[], myInitialAssignment *myInitAssign[], myAlgebraicEquations **myAlgEq, timeVariantAssignments **timeVarAssign, double sim_time, double dt, double *time, allocated_memory *mem, copied_AST *cp_AST){
  unsigned int i, j;
  int flag;
  /* num of each objects */
  unsigned int num_of_species = Model_getNumSpecies(m);
  unsigned int num_of_parameters = Model_getNumParameters(m);
  unsigned int num_of_compartments = Model_getNumCompartments(m);
  unsigned int num_of_reactions = Model_getNumReactions(m);
  unsigned int num_of_rules = Model_getNumRules(m);
  unsigned int num_of_events = Model_getNumEvents(m);
  unsigned int num_of_initialAssignments = Model_getNumInitialAssignments(m);
  unsigned int num_of_time_variant_targets;
  unsigned int num_of_algebraic_rules;

  /* Prepare objects */
  ASTNode_t *node;
  Reaction_t *re;
  InitialAssignment_t *initAssign;
  Rule_t *rule;
  ASTNode_t *times_node, *conv_factor_node;
  Event_t *event;
  Delay_t *delay;

  myEventAssignment *myEvAssign;
  myDelay *evDelay;
  myAlgebraicEquations *algEq;
  char **time_variant_target_id = (char **)malloc(sizeof(char *) * MAX_DELAY_REACTION_NUM);

  /* create mySpecies, myParameters, myCompartments */
  create_species(mySp, m);
  create_parameters(myParam, m);
  create_compartments(myComp, m);

  /* determine species locating compartment */
  locate_species_in_compartment(mySp, num_of_species, myComp, num_of_compartments);

  /* create myReaction & mySpeciseReference without equation */
  create_reactions(myRe, mySp, m);

  /* create myInitialAssignments */
  num_of_time_variant_targets = 0;
  for (i = 0; i < num_of_initialAssignments; i++) {
    myInitAssign[i] = myInitialAssignment_create();
    myInitialAssignment_initWithModel(myInitAssign[i], m, i);
    myInitialAssignment_initTarget(myInitAssign[i], mySp, num_of_species, myParam, num_of_parameters, myComp, num_of_compartments);
    myInitialAssignment_initTargetSpeciesReference(myInitAssign[i], myRe, num_of_reactions);

    initAssign = myInitialAssignment_getOrigin(myInitAssign[i]);
    node = (ASTNode_t*)InitialAssignment_getMath(initAssign);
    node = ASTNode_deepCopy(node);
    TRACE(("original math : "));
    check_AST(node, NULL);
    /* alter_tree_structure(m, &node, cp_AST); */
    alter_tree_structure(m, &node, NULL, 0, cp_AST);
    /* unit */
    if(myInitAssign[i]->target_species != NULL){
      if(myInitAssign[i]->target_species->is_amount
          && !myInitAssign[i]->target_species->has_only_substance_units
          && Compartment_getSpatialDimensions(myInitAssign[i]->target_species->locating_compartment->origin) != 0){
        assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myInitAssign[i]->target_species->origin))), 0);
      }else if(myInitAssign[i]->target_species->is_concentration
          && myInitAssign[i]->target_species->has_only_substance_units){
        assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myInitAssign[i]->target_species->origin))), 1);
      }
    }
    /* unit */
    TRACE(("altered math : "));
    check_AST(node, NULL);
    if(include_time(node, 0)){
      time_variant_target_id[num_of_time_variant_targets++] = (char*)InitialAssignment_getSymbol(initAssign);
    }
    myInitAssign[i]->eq->math_length = get_equation(m, myInitAssign[i]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, NULL, time_variant_target_id, 0, NULL, mem);
    /* ASTNode_free(node); */
    TRACE(("math\n"));
    check_math(myInitAssign[i]->eq);
  }

  /* find time variant target of assignment rule */
  *timeVarAssign = (timeVariantAssignments*)malloc(sizeof(timeVariantAssignments));
  (*timeVarAssign)->num_of_time_variant_assignments = 0;
  for(i=0; i<num_of_rules; i++){
    rule = Model_getRule(m, i);
    node = (ASTNode_t*)Rule_getMath(rule);
    if(Rule_isAssignment(rule) && include_time(node, 0)){
      node = ASTNode_deepCopy(node);
      (*timeVarAssign)->target_id[(*timeVarAssign)->num_of_time_variant_assignments] = (char*)Rule_getVariable(rule);
      TRACE(("original math : "));
      check_AST(node, NULL);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      /* unit */
      for(j=0; j<num_of_species; j++){
        if(strcmp(Species_getId(mySp[j]->origin), Rule_getVariable(rule)) == 0){
          if(mySp[j]->is_amount
              && !mySp[j]->has_only_substance_units
              && Compartment_getSpatialDimensions(mySp[j]->locating_compartment->origin) != 0){
            assignment_alter_tree_structure(&node, (char*)Compartment_getId(mySp[j]->locating_compartment->origin), 0);
          }else if(mySp[j]->is_concentration
              && mySp[j]->has_only_substance_units){
            assignment_alter_tree_structure(&node, (char*)Compartment_getId(mySp[j]->locating_compartment->origin), 1);
          }
          break;
        }
      }
      /* unit */
      TRACE(("altered math : "));
      check_AST(node, NULL);
      (*timeVarAssign)->eq[(*timeVarAssign)->num_of_time_variant_assignments] = equation_create();
      (*timeVarAssign)->eq[(*timeVarAssign)->num_of_time_variant_assignments]->math_length = get_equation(m, (*timeVarAssign)->eq[(*timeVarAssign)->num_of_time_variant_assignments], mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, NULL, time_variant_target_id, 0, NULL, mem);
      ((*timeVarAssign)->num_of_time_variant_assignments)++;
    }
  }

  /* create myReaction & mySpeciseReference equation */
  for(i=0; i<num_of_reactions; i++){
    re = (Reaction_t*)ListOf_get(Model_getListOfReactions(m), i);
    node = (ASTNode_t*)KineticLaw_getMath(Reaction_getKineticLaw(re));
    node = ASTNode_deepCopy(node);
    TRACE(("original math of %s: ", Reaction_getId(myRe[i]->origin)));
    check_AST(node, NULL);
    alter_tree_structure(m, &node, NULL, 0, cp_AST);
    set_local_para_as_value(node, Reaction_getKineticLaw(re));
    TRACE(("alterated math of %s : ", Reaction_getId(myRe[i]->origin)));
    check_AST(node, NULL);
    myRe[i]->eq->math_length = get_equation(m, myRe[i]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
    /* ASTNode_free(node); */
    TRACE(("math of %s\n", Reaction_getId(myRe[i]->origin)));
    check_math(myRe[i]->eq);
    /* products start */
    for(j=0; j<myRe[i]->num_of_products; j++){
      node = (ASTNode_t*)StoichiometryMath_getMath(SpeciesReference_getStoichiometryMath(myRe[i]->products[j]->origin));
      if(node != NULL){/* l2v4 */
        node = ASTNode_deepCopy(node);
        myRe[i]->products[j]->value = 0;
        myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
      }else if(SpeciesReference_isSetStoichiometry(myRe[i]->products[j]->origin) && !SpeciesReference_isSetId(myRe[i]->products[j]->origin)){
        node = ASTNode_createWithType(AST_REAL);
        ASTNode_setReal(node, SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin));
        myRe[i]->products[j]->value = SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin);
        myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
      }else if(SpeciesReference_isSetStoichiometry(myRe[i]->products[j]->origin) && SpeciesReference_isSetId(myRe[i]->products[j]->origin)){
        node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(node, SpeciesReference_getId(myRe[i]->products[j]->origin));
        myRe[i]->products[j]->value = SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin);
        myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
      }
      if(node == NULL){
        if(my_isnan(SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin)) && SpeciesReference_isSetId(myRe[i]->products[j]->origin)){/* l3v1 */
          node = ASTNode_createWithType(AST_NAME);
          ASTNode_setName(node, SpeciesReference_getId(myRe[i]->products[j]->origin));
          myRe[i]->products[j]->value = 1.0;
          myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
        }else if(my_isnan(SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin))){/* l3v1 */
          node = ASTNode_createWithType(AST_REAL);
          ASTNode_setReal(node, 1.0);
          myRe[i]->products[j]->value = 1.0;
          myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
        }else{/* l2v4 */
          node = ASTNode_createWithType(AST_REAL);
          ASTNode_setReal(node, SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin));
          myRe[i]->products[j]->value = SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin);
          myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
        }
      }
      TRACE(("original math of %s : ", SpeciesReference_getSpecies(myRe[i]->products[j]->origin)));
      check_AST(node, NULL);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      /* unit */
      if(myRe[i]->products[j]->mySp->is_concentration){
        assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myRe[i]->products[j]->mySp->origin))), 1);
      }
      /* unit */
      /* conversion factor */
      if(Species_isSetConversionFactor(myRe[i]->products[j]->mySp->origin)){
        times_node = ASTNode_createWithType(AST_TIMES);
        conv_factor_node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(conv_factor_node, Species_getConversionFactor(myRe[i]->products[j]->mySp->origin));
        ASTNode_addChild(times_node, conv_factor_node);
        ASTNode_addChild(times_node, node);
        node = times_node;
      }else if(Model_isSetConversionFactor(m)){
        times_node = ASTNode_createWithType(AST_TIMES);
        conv_factor_node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(conv_factor_node, Model_getConversionFactor(m));
        ASTNode_addChild(times_node, conv_factor_node);
        ASTNode_addChild(times_node, node);
        node = times_node;
      }
      /* conversion factor */
      myRe[i]->products[j]->eq->math_length = get_equation(m, myRe[i]->products[j]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
      TRACE(("alterated math of %s : ", SpeciesReference_getSpecies(myRe[i]->products[j]->origin)));
      check_AST(node, NULL);
      /* ASTNode_free(node); */
      TRACE(("math of %s\n", SpeciesReference_getSpecies(myRe[i]->products[j]->origin)));
      check_math(myRe[i]->products[j]->eq);
      myRe[i]->products[j]->k[0] = 0;
      myRe[i]->products[j]->k[1] = 0;
      myRe[i]->products[j]->k[2] = 0;
      myRe[i]->products[j]->k[3] = 0;
      myRe[i]->products[j]->prev_val[0] = myRe[i]->products[j]->value;
      myRe[i]->products[j]->prev_val[1] = myRe[i]->products[j]->value;
      myRe[i]->products[j]->prev_val[2] = myRe[i]->products[j]->value;
      myRe[i]->products[j]->prev_k[0] = 0;
      myRe[i]->products[j]->prev_k[1] = 0;
      myRe[i]->products[j]->prev_k[2] = 0;
    }
    /* products fin */
    /* reactants start */
    for(j=0; j<myRe[i]->num_of_reactants; j++){
      myRe[i]->reactants[j]->depending_rule = NULL;
      node = (ASTNode_t*)StoichiometryMath_getMath(SpeciesReference_getStoichiometryMath(myRe[i]->reactants[j]->origin));
      if(node != NULL){/* l2v4 */
        node = ASTNode_deepCopy(node);
        myRe[i]->reactants[j]->value = 0;
        myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
      }else if(SpeciesReference_isSetStoichiometry(myRe[i]->reactants[j]->origin) && !SpeciesReference_isSetId(myRe[i]->reactants[j]->origin)){
        node = ASTNode_createWithType(AST_REAL);
        ASTNode_setReal(node, SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin));
        myRe[i]->reactants[j]->value = SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin);
        myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
      }else if(SpeciesReference_isSetStoichiometry(myRe[i]->reactants[j]->origin) && SpeciesReference_isSetId(myRe[i]->reactants[j]->origin)){
        node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(node, SpeciesReference_getId(myRe[i]->reactants[j]->origin));
        myRe[i]->reactants[j]->value = SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin);
        myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
      }
      if(node == NULL){
        if(my_isnan(SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin)) && SpeciesReference_isSetId(myRe[i]->reactants[j]->origin)){/* l3v1 */
          node = ASTNode_createWithType(AST_NAME);
          ASTNode_setName(node, SpeciesReference_getId(myRe[i]->reactants[j]->origin));
          myRe[i]->reactants[j]->value = 1.0;
          myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
        }else if(my_isnan(SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin))){/* l3v1 */
          node = ASTNode_createWithType(AST_REAL);
          ASTNode_setReal(node, 1.0);
          myRe[i]->reactants[j]->value = 1.0;
          myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
        }else{/* l2v4 */
          node = ASTNode_createWithType(AST_REAL);
          ASTNode_setReal(node, SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin));
          myRe[i]->reactants[j]->value = SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin);
          myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
        }
      }
      TRACE(("original math of %s : ", SpeciesReference_getSpecies(myRe[i]->reactants[j]->origin)));
      check_AST(node, NULL);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      /* unit */
      if(myRe[i]->reactants[j]->mySp->is_concentration){
        assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myRe[i]->reactants[j]->mySp->origin))), 1);
      }
      /* unit */
      /* conversion factor */
      if(Species_isSetConversionFactor(myRe[i]->reactants[j]->mySp->origin)){
        times_node = ASTNode_createWithType(AST_TIMES);
        conv_factor_node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(conv_factor_node, Species_getConversionFactor(myRe[i]->reactants[j]->mySp->origin));
        ASTNode_addChild(times_node, conv_factor_node);
        ASTNode_addChild(times_node, node);
        node = times_node;
      }else if(Model_isSetConversionFactor(m)){
        times_node = ASTNode_createWithType(AST_TIMES);
        conv_factor_node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(conv_factor_node, Model_getConversionFactor(m));
        ASTNode_addChild(times_node, conv_factor_node);
        ASTNode_addChild(times_node, node);
        node = times_node;
      }
      /* conversion factor */
      myRe[i]->reactants[j]->eq->math_length = get_equation(m, myRe[i]->reactants[j]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
      TRACE(("altered math of %s : ", SpeciesReference_getSpecies(myRe[i]->reactants[j]->origin)));
      check_AST(node, NULL);
      /* ASTNode_free(node); */
      TRACE(("math of %s\n", SpeciesReference_getSpecies(myRe[i]->reactants[j]->origin)));
      check_math(myRe[i]->reactants[j]->eq);
      myRe[i]->reactants[j]->k[0] = 0;      
      myRe[i]->reactants[j]->k[1] = 0;      
      myRe[i]->reactants[j]->k[2] = 0;      
      myRe[i]->reactants[j]->k[3] = 0;      
      myRe[i]->reactants[j]->prev_val[0] = myRe[i]->reactants[j]->value;
      myRe[i]->reactants[j]->prev_val[1] = myRe[i]->reactants[j]->value;
      myRe[i]->reactants[j]->prev_val[2] = myRe[i]->reactants[j]->value;
      myRe[i]->reactants[j]->prev_k[0] = 0;
      myRe[i]->reactants[j]->prev_k[1] = 0;
      myRe[i]->reactants[j]->prev_k[2] = 0;
      /* reactants fin */
    }
  }

  /* prepare reversible fast reaction */
  prepare_reversible_fast_reaction(m, myRe, mySp, myParam, myComp, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem, cp_AST);

  /* create myRules */
  for (i = 0; i < num_of_rules; i++) {
    myRu[i] = myRule_create();
    myRule_initWithModel(myRu[i], m, i);
    rule = myRule_getOrigin(myRu[i]);
    if (Rule_isRate(rule) || Rule_isAssignment(rule)) {
      myRule_initTarget(myRu[i], mySp, num_of_species, myParam, num_of_parameters, myComp, num_of_compartments, myRe, num_of_reactions);
      node = (ASTNode_t*)Rule_getMath(rule);
      node = ASTNode_deepCopy(node);
      TRACE(("original math : "));
      check_AST(node, NULL);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      /* unit */
      if(myRu[i]->target_species != NULL){
        if(myRu[i]->target_species->is_amount
            && !myRu[i]->target_species->has_only_substance_units
            && Compartment_getSpatialDimensions(myRu[i]->target_species->locating_compartment->origin) != 0){
          assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myRu[i]->target_species->origin))), 0);
        }else if(myRu[i]->target_species->is_concentration
            && myRu[i]->target_species->has_only_substance_units){
          assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myRu[i]->target_species->origin))), 1);
        }
      }
      /* unit */
      TRACE(("altered math : "));
      check_AST(node, NULL);
      myRu[i]->eq->math_length = get_equation(m, myRu[i]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
      /* ASTNode_free(node); */
      TRACE(("math\n"));
      check_math(myRule_getEquation(myRu[i]));
    }
  }

  /* create myEvents & myEventAssignments */
  for (i = 0; i < num_of_events; i++) {
    myEv[i] = myEvent_create();
    myEvent_initWithModel(myEv[i], m, i);
    event = myEvent_getOrigin(myEv[i]);
    node = (ASTNode_t*)Trigger_getMath(Event_getTrigger(event));
    node = ASTNode_deepCopy(node);
    TRACE(("original math of %s : ", Event_getId(myEv[i]->origin)));
    check_AST(node, NULL);
    /* alter_tree_structure(m, &node, cp_AST); */
    alter_tree_structure(m, &node, NULL, 0, cp_AST);
    myEv[i]->eq->math_length = get_equation(m, myEv[i]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
    TRACE(("altered math of %s : ", Event_getId(myEv[i]->origin)));
    check_AST(node, NULL);
    /* ASTNode_free(node); */
    TRACE(("math of %s\n", Event_getId(myEv[i]->origin)));
    check_math(myEv[i]->eq);
    myEv[i]->assignments = (myEventAssignment**)malloc(sizeof(myEventAssignment*)*Event_getNumEventAssignments(event));
    for (j = 0; j < Event_getNumEventAssignments(event); j++) {
      myEvAssign = myEventAssignment_create();
      myEventAssignment_initWithEvent(myEvAssign, event, j);
      myEventAssignment_initTarget(myEvAssign, mySp, num_of_species, myParam, num_of_parameters, myComp, num_of_compartments, myRe, num_of_reactions);
      myEv[i]->assignments[j] = myEvAssign;
      node = (ASTNode_t*)EventAssignment_getMath(myEv[i]->assignments[j]->origin);
      node = ASTNode_deepCopy(node);
      TRACE(("original math : "));
      check_AST(node, NULL);
      if(Event_getDelay(myEv[i]->origin) != NULL && Event_getUseValuesFromTriggerTime(myEv[i]->origin)){
        pre_ev_alter_tree_structure(&node, NULL, 0, (ASTNode_t*)Delay_getMath(Event_getDelay(myEv[i]->origin)));
        TRACE(("after pre ev alter math : "));
        check_AST(node, NULL);
        ev_alter_tree_structure(m, &node, NULL, 0, cp_AST);
        TRACE(("after ev alter math : "));
        check_AST(node, NULL);
        post_ev_alter_tree_structure(m, &node, NULL, 0);
        TRACE(("after post ev alter math : "));
        check_AST(node, NULL);
      }else{
        alter_tree_structure(m, &node, NULL, 0, cp_AST);
      }
      /* unit */
      if(myEv[i]->assignments[j]->target_species != NULL){
        if(myEv[i]->assignments[j]->target_species->is_amount
            && !myEv[i]->assignments[j]->target_species->has_only_substance_units
            && Compartment_getSpatialDimensions(myEv[i]->assignments[j]->target_species->locating_compartment->origin) != 0){
          assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myEv[i]->assignments[j]->target_species->origin))), 0);
        }else if(myEv[i]->assignments[j]->target_species->is_concentration
            && myEv[i]->assignments[j]->target_species->has_only_substance_units){
          assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myEv[i]->assignments[j]->target_species->origin))), 1);
        }
      }
      /* unit */
      TRACE(("altered math : "));
      check_AST(node, NULL);
      myEv[i]->assignments[j]->eq->math_length = get_equation(m, myEv[i]->assignments[j]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
      /* ASTNode_free(node); */
      TRACE(("math\n"));
      check_math(myEv[i]->assignments[j]->eq);
    }
    //myEv[i]->event_delay = NULL;
    delay = Event_getDelay(event);
    if (delay != NULL) {
      evDelay = myDelay_create();
      myDelay_initWithOrigin(evDelay, delay);
      myEv[i]->event_delay = evDelay;
      /*
      myEv[i]->event_delay = (myDelay*)malloc(sizeof(myDelay));
      myEv[i]->event_delay->origin = Event_getDelay(event);
      */
      node = (ASTNode_t *)Delay_getMath(delay);
      node = ASTNode_deepCopy(node);
      TRACE(("original math : "));
      check_AST(node, NULL);
      /* alter_tree_structure(m, &node, cp_AST); */
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      //myEv[i]->event_delay->eq = equation_create();
      myEv[i]->event_delay->eq->math_length = get_equation(m, myEv[i]->event_delay->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
      TRACE(("altered math : "));
      check_AST(node, NULL);
      /* ASTNode_free(node); */
      TRACE(("math\n"));
      check_math(myEv[i]->event_delay->eq);
      myEv[i]->firing_times = (double*)malloc(sizeof(double)*(int)(sim_time/dt));
      for(j=0; j<(unsigned int)(sim_time/dt); j++){
        myEv[i]->firing_times[j] = sim_time+1;
      }
      myEv[i]->num_of_delayed_events_que = 0;
      myEv[i]->next_firing_index = 0;
    }
    //myEv[i]->priority_eq = NULL;
    if(Event_isSetPriority(myEv[i]->origin)){
      myEv[i]->priority_eq = equation_create();
      node = (ASTNode_t*)Priority_getMath(Event_getPriority(myEv[i]->origin));
      node = ASTNode_deepCopy(node);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      myEv[i]->priority_eq->math_length = get_equation(m, myEv[i]->priority_eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
    }
  }

  /* prepare Algebraic */
  num_of_algebraic_rules = 0;
  for(i=0; i<num_of_rules; i++){
    if(myRu[i]->is_algebraic){
      num_of_algebraic_rules++;
    }
  }
  if(num_of_algebraic_rules != 0){
    *myAlgEq = (myAlgebraicEquations*)malloc(sizeof(myAlgebraicEquations));
    algEq = *myAlgEq;
    algEq->num_of_algebraic_rules = num_of_algebraic_rules;
    algEq->coefficient_matrix = NULL;
    algEq->constant_vector = NULL;
    algEq->coefficient = NULL;
    algEq->constant = NULL;
    algEq->target_species = NULL;
    algEq->target_parameter = NULL;
    algEq->target_compartment = NULL;
    for(i=0; i<MAX_ALGEBRAIC_VARIABLES; i++){
      algEq->alg_target_species[i] = NULL;
      algEq->alg_target_parameter[i] = NULL;
      algEq->alg_target_compartment[i] = NULL;
    }
    algEq->num_of_alg_target_sp = 0;
    algEq->num_of_alg_target_param = 0;
    algEq->num_of_alg_target_comp = 0;
    algEq->num_of_algebraic_variables = 0;
    if(algEq->num_of_algebraic_rules > 1){
      algEq->coefficient_matrix = (equation***)malloc(sizeof(equation**)*algEq->num_of_algebraic_rules);
      for(i=0; i<algEq->num_of_algebraic_rules; i++){
        algEq->coefficient_matrix[i] = (equation**)malloc(sizeof(equation*)*algEq->num_of_algebraic_rules);
        for(j=0; j<algEq->num_of_algebraic_rules; j++){
          algEq->coefficient_matrix[i][j] = equation_create();
          algEq->coefficient_matrix[i][j]->math_length = 0;
        }
      }
      algEq->constant_vector = (equation**)malloc(sizeof(equation*)*algEq->num_of_algebraic_rules);
      for(i=0; i<algEq->num_of_algebraic_rules; i++){
        algEq->constant_vector[i] = equation_create();
      }
    }else{
      algEq->coefficient = equation_create();
      algEq->constant = equation_create();
    }
    TRACE(("prepare algebraic start\n"));
    prepare_algebraic(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, algEq, sim_time, dt, time, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem, cp_AST);
    TRACE(("prepare algebraic finish\n"));
    if(algEq->num_of_algebraic_rules > 1){
      for(i=0; i<algEq->num_of_algebraic_rules; i++){
        for(j=0; j<algEq->num_of_algebraic_rules; j++){
          TRACE(("math of coefficient matrix[%d][%d] is\n", i, j));
          check_math(algEq->coefficient_matrix[i][j]);
        }
      }
      for(i=0; i<algEq->num_of_algebraic_rules; i++){
        TRACE(("constant vector[%d] is\n", i));
        check_math(algEq->constant_vector[i]);
      }
    }else{
      TRACE(("math of coefficient is\n"));
      check_math(algEq->coefficient);
      TRACE(("math of constant is\n"));
      check_math(algEq->constant);
    }
    if(algEq->num_of_algebraic_rules > 1){
      flag = 0;
      for(i=0; i<num_of_species; i++){
        for(j=0; j<algEq->num_of_algebraic_variables; j++){
          if(strcmp(algEq->variables_id[j], Species_getId(mySp[i]->origin)) == 0){
            algEq->alg_target_species[algEq->num_of_alg_target_sp] = (myAlgTargetSp*)malloc(sizeof(myAlgTargetSp));
            algEq->alg_target_species[algEq->num_of_alg_target_sp]->target_species = mySp[i];
            algEq->alg_target_species[algEq->num_of_alg_target_sp]->order = j;
            algEq->num_of_alg_target_sp++;
            flag = 1;
            break;
          }
        }
      }
      if(!flag){
        for(i=0; i<num_of_parameters; i++){
          for(j=0; j<algEq->num_of_algebraic_variables; j++){
            if(strcmp(algEq->variables_id[j], Parameter_getId(myParam[i]->origin)) == 0){
              algEq->alg_target_parameter[algEq->num_of_alg_target_param] = (myAlgTargetParam*)malloc(sizeof(myAlgTargetParam));
              algEq->alg_target_parameter[algEq->num_of_alg_target_param]->target_parameter = myParam[i];
              algEq->alg_target_parameter[algEq->num_of_alg_target_param]->order = j;
              algEq->num_of_alg_target_param++;
              flag = 1;
              break;
            }
          }
        }
      }
      if(!flag){
        for(i=0; i<num_of_compartments; i++){
          for(j=0; j<algEq->num_of_algebraic_variables; j++){
            if(strcmp(algEq->variables_id[j], Compartment_getId(myComp[i]->origin)) == 0){
              algEq->alg_target_compartment[algEq->num_of_alg_target_comp] = (myAlgTargetComp*)malloc(sizeof(myAlgTargetComp));
              algEq->alg_target_compartment[algEq->num_of_alg_target_comp]->target_compartment = myComp[i];
              algEq->alg_target_compartment[algEq->num_of_alg_target_comp]->order = j;
              algEq->num_of_alg_target_comp++;
              flag = 1;
              break;
            }
          }
        }
      }
    }else{
      flag = 0;
      for(i=0; i<num_of_species; i++){
        if(strcmp(algEq->variables_id[0], Species_getId(mySp[i]->origin)) == 0){
          algEq->target_species = mySp[i];
          flag = 1;
          break;
        }
      }
      if(!flag){
        for(i=0; i<num_of_parameters; i++){
          if(strcmp(algEq->variables_id[0], Parameter_getId(myParam[i]->origin)) == 0){
            algEq->target_parameter = myParam[i];
            flag = 1;
            break;
          }
        }
      }
      if(!flag){
        for(i=0; i<num_of_compartments; i++){
          if(strcmp(algEq->variables_id[0], Compartment_getId(myComp[i]->origin)) == 0){
            algEq->target_compartment = myComp[i];
            flag = 1;
            break;
          }
        }
      }
    }
  }
}

/* create my SBML obejects for efficient simulations (for bifurcation analysis, instructions on local parameter in Kinetic Law are a little changed)*/
void create_mySBML_objects_forBA(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], myReaction *myRe[], myRule *myRu[], myEvent *myEv[], myInitialAssignment *myInitAssign[], myAlgebraicEquations **myAlgEq, timeVariantAssignments **timeVarAssign, double sim_time, double dt, double *time, allocated_memory *mem, copied_AST *cp_AST, char* bif_param_id, double bif_param_value){
  unsigned int i, j;
  int flag;
  /* num of each objects */
  unsigned int num_of_species = Model_getNumSpecies(m);
  unsigned int num_of_parameters = Model_getNumParameters(m);
  unsigned int num_of_compartments = Model_getNumCompartments(m);
  unsigned int num_of_reactions = Model_getNumReactions(m);
  unsigned int num_of_rules = Model_getNumRules(m);
  unsigned int num_of_events = Model_getNumEvents(m);
  unsigned int num_of_initialAssignments = Model_getNumInitialAssignments(m);
  unsigned int num_of_time_variant_targets;
  unsigned int num_of_algebraic_rules;

  /* Prepare objects */
  ASTNode_t *node;
  Reaction_t *re;
  InitialAssignment_t *initAssign;
  Rule_t *rule;
  ASTNode_t *times_node, *conv_factor_node;
  Event_t *event;

  myEventAssignment *myEvAssign;
  myAlgebraicEquations *algEq;
  char **time_variant_target_id = (char **)malloc(sizeof(char *) * MAX_DELAY_REACTION_NUM);

  /* create mySpecies, myParameters, myCompartments */
  create_species(mySp, m);
  create_parameters(myParam, m);
  create_compartments(myComp, m);

  /* determine species locating compartment */
  locate_species_in_compartment(mySp, num_of_species, myComp, num_of_compartments);

  /* create myReaction & mySpeciseReference without equation */
  create_reactions(myRe, mySp, m);

  /* create myInitialAssignments */
  num_of_time_variant_targets = 0;
  for (i = 0; i < num_of_initialAssignments; i++) {
    myInitAssign[i] = myInitialAssignment_create();
    myInitialAssignment_initWithModel(myInitAssign[i], m, i);
    myInitialAssignment_initTarget(myInitAssign[i], mySp, num_of_species, myParam, num_of_parameters, myComp, num_of_compartments);
    myInitialAssignment_initTargetSpeciesReference(myInitAssign[i], myRe, num_of_reactions);

    initAssign = myInitialAssignment_getOrigin(myInitAssign[i]);
    node = (ASTNode_t*)InitialAssignment_getMath(initAssign);
    node = ASTNode_deepCopy(node);
    TRACE(("original math : "));
    check_AST(node, NULL);
    /* alter_tree_structure(m, &node, cp_AST); */
    alter_tree_structure(m, &node, NULL, 0, cp_AST);
    /* unit */
    if(myInitAssign[i]->target_species != NULL){
      if(myInitAssign[i]->target_species->is_amount
          && !myInitAssign[i]->target_species->has_only_substance_units
          && Compartment_getSpatialDimensions(myInitAssign[i]->target_species->locating_compartment->origin) != 0){
        assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myInitAssign[i]->target_species->origin))), 0);
      }else if(myInitAssign[i]->target_species->is_concentration
          && myInitAssign[i]->target_species->has_only_substance_units){
        assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myInitAssign[i]->target_species->origin))), 1);
      }
    }
    /* unit */
    TRACE(("altered math : "));
    check_AST(node, NULL);
    if(include_time(node, 0)){
      time_variant_target_id[num_of_time_variant_targets++] = (char*)InitialAssignment_getSymbol(initAssign);
    }
    myInitAssign[i]->eq->math_length = get_equation(m, myInitAssign[i]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, NULL, time_variant_target_id, 0, NULL, mem);
    /* ASTNode_free(node); */
    TRACE(("math\n"));
    check_math(myInitAssign[i]->eq);
  }

  /* find time variant target of assignment rule */
  *timeVarAssign = (timeVariantAssignments*)malloc(sizeof(timeVariantAssignments));
  (*timeVarAssign)->num_of_time_variant_assignments = 0;
  for(i=0; i<num_of_rules; i++){
    rule = Model_getRule(m, i);
    node = (ASTNode_t*)Rule_getMath(rule);
    if(Rule_isAssignment(rule) && include_time(node, 0)){
      node = ASTNode_deepCopy(node);
      (*timeVarAssign)->target_id[(*timeVarAssign)->num_of_time_variant_assignments] = (char*)Rule_getVariable(rule);
      TRACE(("original math : "));
      check_AST(node, NULL);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      /* unit */
      for(j=0; j<num_of_species; j++){
        if(strcmp(Species_getId(mySp[j]->origin), Rule_getVariable(rule)) == 0){
          if(mySp[j]->is_amount
              && !mySp[j]->has_only_substance_units
              && Compartment_getSpatialDimensions(mySp[j]->locating_compartment->origin) != 0){
            assignment_alter_tree_structure(&node, (char*)Compartment_getId(mySp[j]->locating_compartment->origin), 0);
          }else if(mySp[j]->is_concentration
              && mySp[j]->has_only_substance_units){
            assignment_alter_tree_structure(&node, (char*)Compartment_getId(mySp[j]->locating_compartment->origin), 1);
          }
          break;
        }
      }
      /* unit */
      TRACE(("altered math : "));
      check_AST(node, NULL);
      (*timeVarAssign)->eq[(*timeVarAssign)->num_of_time_variant_assignments] = equation_create();
      (*timeVarAssign)->eq[(*timeVarAssign)->num_of_time_variant_assignments]->math_length = get_equation(m, (*timeVarAssign)->eq[(*timeVarAssign)->num_of_time_variant_assignments], mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, NULL, time_variant_target_id, 0, NULL, mem);
      ((*timeVarAssign)->num_of_time_variant_assignments)++;
    }
  }

  /* create myReaction & mySpeciseReference equation */
  for(i=0; i<num_of_reactions; i++){
    re = (Reaction_t*)ListOf_get(Model_getListOfReactions(m), i);
    node = (ASTNode_t*)KineticLaw_getMath(Reaction_getKineticLaw(re));
    node = ASTNode_deepCopy(node);
    TRACE(("original math of %s: ", Reaction_getId(myRe[i]->origin)));
    check_AST(node, NULL);
    alter_tree_structure(m, &node, NULL, 0, cp_AST);
	set_local_para_as_value_forBA(node, Reaction_getKineticLaw(re), bif_param_id, bif_param_value);
	TRACE(("alterated math of %s : ", Reaction_getId(myRe[i]->origin)));
    check_AST(node, NULL);
    myRe[i]->eq->math_length = get_equation(m, myRe[i]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
    /* ASTNode_free(node); */
    TRACE(("math of %s\n", Reaction_getId(myRe[i]->origin)));
    check_math(myRe[i]->eq);
    /* products start */
    for(j=0; j<myRe[i]->num_of_products; j++){
      node = (ASTNode_t*)StoichiometryMath_getMath(SpeciesReference_getStoichiometryMath(myRe[i]->products[j]->origin));
      if(node != NULL){/* l2v4 */
        node = ASTNode_deepCopy(node);
        myRe[i]->products[j]->value = 0;
        myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
      }else if(SpeciesReference_isSetStoichiometry(myRe[i]->products[j]->origin) && !SpeciesReference_isSetId(myRe[i]->products[j]->origin)){
        node = ASTNode_createWithType(AST_REAL);
        ASTNode_setReal(node, SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin));
        myRe[i]->products[j]->value = SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin);
        myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
      }else if(SpeciesReference_isSetStoichiometry(myRe[i]->products[j]->origin) && SpeciesReference_isSetId(myRe[i]->products[j]->origin)){
        node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(node, SpeciesReference_getId(myRe[i]->products[j]->origin));
        myRe[i]->products[j]->value = SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin);
        myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
      }
      if(node == NULL){
        if(my_isnan(SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin)) && SpeciesReference_isSetId(myRe[i]->products[j]->origin)){/* l3v1 */
          node = ASTNode_createWithType(AST_NAME);
          ASTNode_setName(node, SpeciesReference_getId(myRe[i]->products[j]->origin));
          myRe[i]->products[j]->value = 1.0;
          myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
        }else if(my_isnan(SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin))){/* l3v1 */
          node = ASTNode_createWithType(AST_REAL);
          ASTNode_setReal(node, 1.0);
          myRe[i]->products[j]->value = 1.0;
          myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
        }else{/* l2v4 */
          node = ASTNode_createWithType(AST_REAL);
          ASTNode_setReal(node, SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin));
          myRe[i]->products[j]->value = SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin);
          myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
        }
      }
      TRACE(("original math of %s : ", SpeciesReference_getSpecies(myRe[i]->products[j]->origin)));
      check_AST(node, NULL);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      /* unit */
      if(myRe[i]->products[j]->mySp->is_concentration){
        assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myRe[i]->products[j]->mySp->origin))), 1);
      }
      /* unit */
      /* conversion factor */
      if(Species_isSetConversionFactor(myRe[i]->products[j]->mySp->origin)){
        times_node = ASTNode_createWithType(AST_TIMES);
        conv_factor_node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(conv_factor_node, Species_getConversionFactor(myRe[i]->products[j]->mySp->origin));
        ASTNode_addChild(times_node, conv_factor_node);
        ASTNode_addChild(times_node, node);
        node = times_node;
      }else if(Model_isSetConversionFactor(m)){
        times_node = ASTNode_createWithType(AST_TIMES);
        conv_factor_node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(conv_factor_node, Model_getConversionFactor(m));
        ASTNode_addChild(times_node, conv_factor_node);
        ASTNode_addChild(times_node, node);
        node = times_node;
      }
      /* conversion factor */
      myRe[i]->products[j]->eq->math_length = get_equation(m, myRe[i]->products[j]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
      TRACE(("alterated math of %s : ", SpeciesReference_getSpecies(myRe[i]->products[j]->origin)));
      check_AST(node, NULL);
      /* ASTNode_free(node); */
      TRACE(("math of %s\n", SpeciesReference_getSpecies(myRe[i]->products[j]->origin)));
      check_math(myRe[i]->products[j]->eq);
      myRe[i]->products[j]->k[0] = 0;
      myRe[i]->products[j]->k[1] = 0;
      myRe[i]->products[j]->k[2] = 0;
      myRe[i]->products[j]->k[3] = 0;

      myRe[i]->products[j]->prev_val[0] = myRe[i]->products[j]->value;
      myRe[i]->products[j]->prev_val[1] = myRe[i]->products[j]->value;
      myRe[i]->products[j]->prev_val[2] = myRe[i]->products[j]->value;
      myRe[i]->products[j]->prev_k[0] = 0;
      myRe[i]->products[j]->prev_k[1] = 0;
      myRe[i]->products[j]->prev_k[2] = 0;
    }
    /* products fin */
    /* reactants start */
    for(j=0; j<myRe[i]->num_of_reactants; j++){
      myRe[i]->reactants[j]->depending_rule = NULL;
      node = (ASTNode_t*)StoichiometryMath_getMath(SpeciesReference_getStoichiometryMath(myRe[i]->reactants[j]->origin));
      if(node != NULL){/* l2v4 */
        node = ASTNode_deepCopy(node);
        myRe[i]->reactants[j]->value = 0;
        myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
      }else if(SpeciesReference_isSetStoichiometry(myRe[i]->reactants[j]->origin) && !SpeciesReference_isSetId(myRe[i]->reactants[j]->origin)){
        node = ASTNode_createWithType(AST_REAL);
        ASTNode_setReal(node, SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin));
        myRe[i]->reactants[j]->value = SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin);
        myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
      }else if(SpeciesReference_isSetStoichiometry(myRe[i]->reactants[j]->origin) && SpeciesReference_isSetId(myRe[i]->reactants[j]->origin)){
        node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(node, SpeciesReference_getId(myRe[i]->reactants[j]->origin));
        myRe[i]->reactants[j]->value = SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin);
        myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
      }
      if(node == NULL){
        if(my_isnan(SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin)) && SpeciesReference_isSetId(myRe[i]->reactants[j]->origin)){/* l3v1 */
          node = ASTNode_createWithType(AST_NAME);
          ASTNode_setName(node, SpeciesReference_getId(myRe[i]->reactants[j]->origin));
          myRe[i]->reactants[j]->value = 1.0;
          myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
        }else if(my_isnan(SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin))){/* l3v1 */
          node = ASTNode_createWithType(AST_REAL);
          ASTNode_setReal(node, 1.0);
          myRe[i]->reactants[j]->value = 1.0;
          myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
        }else{/* l2v4 */
          node = ASTNode_createWithType(AST_REAL);
          ASTNode_setReal(node, SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin));
          myRe[i]->reactants[j]->value = SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin);
          myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
        }
      }
      TRACE(("original math of %s : ", SpeciesReference_getSpecies(myRe[i]->reactants[j]->origin)));
      check_AST(node, NULL);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      /* unit */
      if(myRe[i]->reactants[j]->mySp->is_concentration){
        assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myRe[i]->reactants[j]->mySp->origin))), 1);
      }
      /* unit */
      /* conversion factor */
      if(Species_isSetConversionFactor(myRe[i]->reactants[j]->mySp->origin)){
        times_node = ASTNode_createWithType(AST_TIMES);
        conv_factor_node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(conv_factor_node, Species_getConversionFactor(myRe[i]->reactants[j]->mySp->origin));
        ASTNode_addChild(times_node, conv_factor_node);
        ASTNode_addChild(times_node, node);
        node = times_node;
      }else if(Model_isSetConversionFactor(m)){
        times_node = ASTNode_createWithType(AST_TIMES);
        conv_factor_node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(conv_factor_node, Model_getConversionFactor(m));
        ASTNode_addChild(times_node, conv_factor_node);
        ASTNode_addChild(times_node, node);
        node = times_node;
      }
      /* conversion factor */
      myRe[i]->reactants[j]->eq->math_length = get_equation(m, myRe[i]->reactants[j]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
      TRACE(("altered math of %s : ", SpeciesReference_getSpecies(myRe[i]->reactants[j]->origin)));
      check_AST(node, NULL);
      /* ASTNode_free(node); */
      TRACE(("math of %s\n", SpeciesReference_getSpecies(myRe[i]->reactants[j]->origin)));
      check_math(myRe[i]->reactants[j]->eq);
      myRe[i]->reactants[j]->k[0] = 0;      
      myRe[i]->reactants[j]->k[1] = 0;      
      myRe[i]->reactants[j]->k[2] = 0;      
      myRe[i]->reactants[j]->k[3] = 0;      

      myRe[i]->reactants[j]->prev_val[0] = myRe[i]->reactants[j]->value;
      myRe[i]->reactants[j]->prev_val[1] = myRe[i]->reactants[j]->value;
      myRe[i]->reactants[j]->prev_val[2] = myRe[i]->reactants[j]->value;
      myRe[i]->reactants[j]->prev_k[0] = 0;
      myRe[i]->reactants[j]->prev_k[1] = 0;
      myRe[i]->reactants[j]->prev_k[2] = 0;
      /* reactants fin */
    }
  }

  /* prepare reversible fast reaction */
  prepare_reversible_fast_reaction(m, myRe, mySp, myParam, myComp, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem, cp_AST);

  /* create myRules */
  for(i=0; i<num_of_rules; i++){
    myRu[i] = myRule_create();
    myRule_initWithModel(myRu[i], m, i);
    rule = myRule_getOrigin(myRu[i]);
    if (Rule_isRate(rule) || Rule_isAssignment(rule)) {
      myRule_initTarget(myRu[i], mySp, num_of_species, myParam, num_of_parameters, myComp, num_of_compartments, myRe, num_of_reactions);
      node = (ASTNode_t*)Rule_getMath(rule);
      node = ASTNode_deepCopy(node);
      TRACE(("original math : "));
      check_AST(node, NULL);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      /* unit */
      if(myRu[i]->target_species != NULL){
        if(myRu[i]->target_species->is_amount
            && !myRu[i]->target_species->has_only_substance_units
            && Compartment_getSpatialDimensions(myRu[i]->target_species->locating_compartment->origin) != 0){
          assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myRu[i]->target_species->origin))), 0);
        }else if(myRu[i]->target_species->is_concentration
            && myRu[i]->target_species->has_only_substance_units){
          assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myRu[i]->target_species->origin))), 1);
        }
      }
      /* unit */
      TRACE(("altered math : "));
      check_AST(node, NULL);
      myRu[i]->eq->math_length = get_equation(m, myRu[i]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
      /* ASTNode_free(node); */
      TRACE(("math\n"));
      check_math(myRule_getEquation(myRu[i]));
    }
  }

  /* create myEvents & myEventAssignments */
  for (i = 0; i < num_of_events; i++) {
    myEv[i] = myEvent_create();
    myEvent_initWithModel(myEv[i], m, i);
    event = myEvent_getOrigin(myEv[i]);
    node = (ASTNode_t*)Trigger_getMath(Event_getTrigger(event));
    node = ASTNode_deepCopy(node);
    TRACE(("original math of %s : ", Event_getId(myEv[i]->origin)));
    check_AST(node, NULL);
    /* alter_tree_structure(m, &node, cp_AST); */
    alter_tree_structure(m, &node, NULL, 0, cp_AST);
    myEv[i]->eq->math_length = get_equation(m, myEv[i]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
    TRACE(("altered math of %s : ", Event_getId(myEv[i]->origin)));
    check_AST(node, NULL);
    /* ASTNode_free(node); */
    TRACE(("math of %s\n", Event_getId(myEv[i]->origin)));
    check_math(myEv[i]->eq);
    myEv[i]->assignments = (myEventAssignment**)malloc(sizeof(myEventAssignment*)*Event_getNumEventAssignments(event));
    for (j = 0; j < Event_getNumEventAssignments(event); j++) {
      myEvAssign = myEventAssignment_create();
      myEventAssignment_initWithEvent(myEvAssign, event, j);
      myEventAssignment_initTarget(myEvAssign, mySp, num_of_species, myParam, num_of_parameters, myComp, num_of_compartments, myRe, num_of_reactions);
      myEv[i]->assignments[j] = myEvAssign;
      node = (ASTNode_t*)EventAssignment_getMath(myEv[i]->assignments[j]->origin);
      node = ASTNode_deepCopy(node);
      TRACE(("original math : "));
      check_AST(node, NULL);
      if(Event_getDelay(myEv[i]->origin) != NULL && Event_getUseValuesFromTriggerTime(myEv[i]->origin)){
        pre_ev_alter_tree_structure(&node, NULL, 0, (ASTNode_t*)Delay_getMath(Event_getDelay(myEv[i]->origin)));
        TRACE(("after pre ev alter math : "));
        check_AST(node, NULL);
        ev_alter_tree_structure(m, &node, NULL, 0, cp_AST);
        TRACE(("after ev alter math : "));
        check_AST(node, NULL);
        post_ev_alter_tree_structure(m, &node, NULL, 0);
        TRACE(("after post ev alter math : "));
        check_AST(node, NULL);
      }else{
        alter_tree_structure(m, &node, NULL, 0, cp_AST);
      }
      /* unit */
      if(myEv[i]->assignments[j]->target_species != NULL){
        if(myEv[i]->assignments[j]->target_species->is_amount
            && !myEv[i]->assignments[j]->target_species->has_only_substance_units
            && Compartment_getSpatialDimensions(myEv[i]->assignments[j]->target_species->locating_compartment->origin) != 0){
          assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myEv[i]->assignments[j]->target_species->origin))), 0);
        }else if(myEv[i]->assignments[j]->target_species->is_concentration
            && myEv[i]->assignments[j]->target_species->has_only_substance_units){
          assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myEv[i]->assignments[j]->target_species->origin))), 1);
        }
      }
      /* unit */
      TRACE(("altered math : "));
      check_AST(node, NULL);
      myEv[i]->assignments[j]->eq->math_length = get_equation(m, myEv[i]->assignments[j]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
      /* ASTNode_free(node); */
      TRACE(("math\n"));
      check_math(myEv[i]->assignments[j]->eq);
    }
    myEv[i]->event_delay = NULL;
    if(Event_getDelay(event) != NULL){
      myEv[i]->event_delay = (myDelay*)malloc(sizeof(myDelay));
      myEv[i]->event_delay->origin = Event_getDelay(event);
      node = (ASTNode_t*)Delay_getMath(myEv[i]->event_delay->origin);
      node = ASTNode_deepCopy(node);
      TRACE(("original math : "));
      check_AST(node, NULL);
      /* alter_tree_structure(m, &node, cp_AST); */
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      myEv[i]->event_delay->eq = equation_create();
      myEv[i]->event_delay->eq->math_length = get_equation(m, myEv[i]->event_delay->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
      TRACE(("altered math : "));
      check_AST(node, NULL);
      /* ASTNode_free(node); */
      TRACE(("math\n"));
      check_math(myEv[i]->event_delay->eq);
      myEv[i]->firing_times = (double*)malloc(sizeof(double)*(int)(sim_time/dt));
      for(j=0; j<(unsigned int)(sim_time/dt); j++){
        myEv[i]->firing_times[j] = sim_time+1;
      }
      myEv[i]->num_of_delayed_events_que = 0;
      myEv[i]->next_firing_index = 0;
    }
    myEv[i]->priority_eq = NULL;
    if(Event_isSetPriority(myEv[i]->origin)){
      myEv[i]->priority_eq = equation_create();
      node = (ASTNode_t*)Priority_getMath(Event_getPriority(myEv[i]->origin));
      node = ASTNode_deepCopy(node);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      myEv[i]->priority_eq->math_length = get_equation(m, myEv[i]->priority_eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
    }
  }

  /* prepare Algebraic */
  num_of_algebraic_rules = 0;
  for(i=0; i<num_of_rules; i++){
    if(myRu[i]->is_algebraic){
      num_of_algebraic_rules++;
    }
  }
  if(num_of_algebraic_rules != 0){
    *myAlgEq = (myAlgebraicEquations*)malloc(sizeof(myAlgebraicEquations));
    algEq = *myAlgEq;
    algEq->num_of_algebraic_rules = num_of_algebraic_rules;
    algEq->coefficient_matrix = NULL;
    algEq->constant_vector = NULL;
    algEq->coefficient = NULL;
    algEq->constant = NULL;
    algEq->target_species = NULL;
    algEq->target_parameter = NULL;
    algEq->target_compartment = NULL;
    for(i=0; i<MAX_ALGEBRAIC_VARIABLES; i++){
      algEq->alg_target_species[i] = NULL;
      algEq->alg_target_parameter[i] = NULL;
      algEq->alg_target_compartment[i] = NULL;
    }
    algEq->num_of_alg_target_sp = 0;
    algEq->num_of_alg_target_param = 0;
    algEq->num_of_alg_target_comp = 0;
    algEq->num_of_algebraic_variables = 0;
    if(algEq->num_of_algebraic_rules > 1){
      algEq->coefficient_matrix = (equation***)malloc(sizeof(equation**)*algEq->num_of_algebraic_rules);
      for(i=0; i<algEq->num_of_algebraic_rules; i++){
        algEq->coefficient_matrix[i] = (equation**)malloc(sizeof(equation*)*algEq->num_of_algebraic_rules);
        for(j=0; j<algEq->num_of_algebraic_rules; j++){
          algEq->coefficient_matrix[i][j] = equation_create();
          algEq->coefficient_matrix[i][j]->math_length = 0;
        }
      }
      algEq->constant_vector = (equation**)malloc(sizeof(equation*)*algEq->num_of_algebraic_rules);
      for(i=0; i<algEq->num_of_algebraic_rules; i++){
        algEq->constant_vector[i] = equation_create();
      }
    }else{
      algEq->coefficient = equation_create();
      algEq->constant = equation_create();
    }
    TRACE(("prepare algebraic start\n"));
    prepare_algebraic(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, algEq, sim_time, dt, time, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem, cp_AST);
    TRACE(("prepare algebraic finish\n"));
    if(algEq->num_of_algebraic_rules > 1){
      for(i=0; i<algEq->num_of_algebraic_rules; i++){
        for(j=0; j<algEq->num_of_algebraic_rules; j++){
          TRACE(("math of coefficient matrix[%d][%d] is\n", i, j));
          check_math(algEq->coefficient_matrix[i][j]);
        }
      }
      for(i=0; i<algEq->num_of_algebraic_rules; i++){
        TRACE(("constant vector[%d] is\n", i));
        check_math(algEq->constant_vector[i]);
      }
    }else{
      TRACE(("math of coefficient is\n"));
      check_math(algEq->coefficient);
      TRACE(("math of constant is\n"));
      check_math(algEq->constant);
    }
    if(algEq->num_of_algebraic_rules > 1){
      flag = 0;
      for(i=0; i<num_of_species; i++){
        for(j=0; j<algEq->num_of_algebraic_variables; j++){
          if(strcmp(algEq->variables_id[j], Species_getId(mySp[i]->origin)) == 0){
            algEq->alg_target_species[algEq->num_of_alg_target_sp] = (myAlgTargetSp*)malloc(sizeof(myAlgTargetSp));
            algEq->alg_target_species[algEq->num_of_alg_target_sp]->target_species = mySp[i];
            algEq->alg_target_species[algEq->num_of_alg_target_sp]->order = j;
            algEq->num_of_alg_target_sp++;
            flag = 1;
            break;
          }
        }
      }
      if(!flag){
        for(i=0; i<num_of_parameters; i++){
          for(j=0; j<algEq->num_of_algebraic_variables; j++){
            if(strcmp(algEq->variables_id[j], Parameter_getId(myParam[i]->origin)) == 0){
              algEq->alg_target_parameter[algEq->num_of_alg_target_param] = (myAlgTargetParam*)malloc(sizeof(myAlgTargetParam));
              algEq->alg_target_parameter[algEq->num_of_alg_target_param]->target_parameter = myParam[i];
              algEq->alg_target_parameter[algEq->num_of_alg_target_param]->order = j;
              algEq->num_of_alg_target_param++;
              flag = 1;
              break;
            }
          }
        }
      }
      if(!flag){
        for(i=0; i<num_of_compartments; i++){
          for(j=0; j<algEq->num_of_algebraic_variables; j++){
            if(strcmp(algEq->variables_id[j], Compartment_getId(myComp[i]->origin)) == 0){
              algEq->alg_target_compartment[algEq->num_of_alg_target_comp] = (myAlgTargetComp*)malloc(sizeof(myAlgTargetComp));
              algEq->alg_target_compartment[algEq->num_of_alg_target_comp]->target_compartment = myComp[i];
              algEq->alg_target_compartment[algEq->num_of_alg_target_comp]->order = j;
              algEq->num_of_alg_target_comp++;
              flag = 1;
              break;
            }
          }
        }
      }
    }else{
      flag = 0;
      for(i=0; i<num_of_species; i++){
        if(strcmp(algEq->variables_id[0], Species_getId(mySp[i]->origin)) == 0){
          algEq->target_species = mySp[i];
          flag = 1;
          break;
        }
      }
      if(!flag){
        for(i=0; i<num_of_parameters; i++){
          if(strcmp(algEq->variables_id[0], Parameter_getId(myParam[i]->origin)) == 0){
            algEq->target_parameter = myParam[i];
            flag = 1;
            break;
          }
        }
      }
      if(!flag){
        for(i=0; i<num_of_compartments; i++){
          if(strcmp(algEq->variables_id[0], Compartment_getId(myComp[i]->origin)) == 0){
            algEq->target_compartment = myComp[i];
            flag = 1;
            break;
          }
        }
      }
    }
  }
}

/* create my SBML obejects for efficient simulations [variable step-size, allocated memory of delay_val is changed (in get_equation), initialize k(for 6 stage)]*/
void create_mySBML_objectsf(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], myReaction *myRe[], myRule *myRu[], myEvent *myEv[], myInitialAssignment *myInitAssign[], myAlgebraicEquations **myAlgEq, timeVariantAssignments **timeVarAssign, double sim_time, double dt, double *time, allocated_memory *mem, copied_AST *cp_AST, int print_interval){
  unsigned int i, j;
  int flag;
  /* num of each objects */
  unsigned int num_of_species = Model_getNumSpecies(m);
  unsigned int num_of_parameters = Model_getNumParameters(m);
  unsigned int num_of_compartments = Model_getNumCompartments(m);
  unsigned int num_of_reactions = Model_getNumReactions(m);
  unsigned int num_of_rules = Model_getNumRules(m);
  unsigned int num_of_events = Model_getNumEvents(m);
  unsigned int num_of_initialAssignments = Model_getNumInitialAssignments(m);
  unsigned int num_of_time_variant_targets;
  unsigned int num_of_algebraic_rules;

  /* Prepare objects */
  ASTNode_t *node;
  Reaction_t *re;
  InitialAssignment_t *initAssign;
  Rule_t *rule;
  ASTNode_t *times_node, *conv_factor_node;
  Event_t *event;

  myEventAssignment *myEvAssign;
  myAlgebraicEquations *algEq;
  char **time_variant_target_id = (char **)malloc(sizeof(char *) * MAX_DELAY_REACTION_NUM);

  /* create mySpecies, myParameters, myCompartments */
  create_species(mySp, m);
  create_parameters(myParam, m);
  create_compartments(myComp, m);

  /* determine species locating compartment */
  locate_species_in_compartment(mySp, num_of_species, myComp, num_of_compartments);

  /* create myReaction & mySpeciseReference without equation */
  create_reactions(myRe, mySp, m);

  /* create myInitialAssignments */
  num_of_time_variant_targets = 0;
  for (i = 0; i < num_of_initialAssignments; i++) {
    myInitAssign[i] = myInitialAssignment_create();
    myInitialAssignment_initWithModel(myInitAssign[i], m, i);
    myInitialAssignment_initTarget(myInitAssign[i], mySp, num_of_species, myParam, num_of_parameters, myComp, num_of_compartments);
    myInitialAssignment_initTargetSpeciesReference(myInitAssign[i], myRe, num_of_reactions);

    initAssign = myInitialAssignment_getOrigin(myInitAssign[i]);
    node = (ASTNode_t*)InitialAssignment_getMath(initAssign);
    node = ASTNode_deepCopy(node);
    TRACE(("original math : "));
    check_AST(node, NULL);
    /* alter_tree_structure(m, &node, cp_AST); */
    alter_tree_structure(m, &node, NULL, 0, cp_AST);
    /* unit */
    if(myInitAssign[i]->target_species != NULL){
      if(myInitAssign[i]->target_species->is_amount
          && !myInitAssign[i]->target_species->has_only_substance_units
          && Compartment_getSpatialDimensions(myInitAssign[i]->target_species->locating_compartment->origin) != 0){
        assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myInitAssign[i]->target_species->origin))), 0);
      }else if(myInitAssign[i]->target_species->is_concentration
          && myInitAssign[i]->target_species->has_only_substance_units){
        assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myInitAssign[i]->target_species->origin))), 1);
      }
    }
    /* unit */
    TRACE(("altered math : "));
    check_AST(node, NULL);
    if(include_time(node, 0)){
      time_variant_target_id[num_of_time_variant_targets++] = (char*)InitialAssignment_getSymbol(initAssign);
    }
    myInitAssign[i]->eq->math_length = get_equationf(m, myInitAssign[i]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, NULL, time_variant_target_id, 0, NULL, mem, print_interval);
    /* ASTNode_free(node); */
    TRACE(("math\n"));
    check_math(myInitAssign[i]->eq);
  }

  /* find time variant target of assignment rule */
  *timeVarAssign = (timeVariantAssignments*)malloc(sizeof(timeVariantAssignments));
  (*timeVarAssign)->num_of_time_variant_assignments = 0;
  for(i=0; i<num_of_rules; i++){
    rule = Model_getRule(m, i);
    node = (ASTNode_t*)Rule_getMath(rule);
    if(Rule_isAssignment(rule) && include_time(node, 0)){
      node = ASTNode_deepCopy(node);
      (*timeVarAssign)->target_id[(*timeVarAssign)->num_of_time_variant_assignments] = (char*)Rule_getVariable(rule);
      TRACE(("original math : "));
      check_AST(node, NULL);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      /* unit */
      for(j=0; j<num_of_species; j++){
        if(strcmp(Species_getId(mySp[j]->origin), Rule_getVariable(rule)) == 0){
          if(mySp[j]->is_amount
              && !mySp[j]->has_only_substance_units
              && Compartment_getSpatialDimensions(mySp[j]->locating_compartment->origin) != 0){
            assignment_alter_tree_structure(&node, (char*)Compartment_getId(mySp[j]->locating_compartment->origin), 0);
          }else if(mySp[j]->is_concentration
              && mySp[j]->has_only_substance_units){
            assignment_alter_tree_structure(&node, (char*)Compartment_getId(mySp[j]->locating_compartment->origin), 1);
          }
          break;
        }
      }
      /* unit */
      TRACE(("altered math : "));
      check_AST(node, NULL);
      (*timeVarAssign)->eq[(*timeVarAssign)->num_of_time_variant_assignments] = equation_create();
      (*timeVarAssign)->eq[(*timeVarAssign)->num_of_time_variant_assignments]->math_length = get_equationf(m, (*timeVarAssign)->eq[(*timeVarAssign)->num_of_time_variant_assignments], mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, NULL, time_variant_target_id, 0, NULL, mem, print_interval);
      ((*timeVarAssign)->num_of_time_variant_assignments)++;
    }
  }

  /* create myReaction & mySpeciseReference equation */
  for(i=0; i<num_of_reactions; i++){
    re = (Reaction_t*)ListOf_get(Model_getListOfReactions(m), i);
    node = (ASTNode_t*)KineticLaw_getMath(Reaction_getKineticLaw(re));
    node = ASTNode_deepCopy(node);
    TRACE(("original math of %s: ", Reaction_getId(myRe[i]->origin)));
    check_AST(node, NULL);
    alter_tree_structure(m, &node, NULL, 0, cp_AST);
	set_local_para_as_value(node, Reaction_getKineticLaw(re));
	TRACE(("alterated math of %s : ", Reaction_getId(myRe[i]->origin)));
    check_AST(node, NULL);
    myRe[i]->eq->math_length = get_equationf(m, myRe[i]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem, print_interval);
    /* ASTNode_free(node); */
    TRACE(("math of %s\n", Reaction_getId(myRe[i]->origin)));
    check_math(myRe[i]->eq);
    /* products start */
    for(j=0; j<myRe[i]->num_of_products; j++){
      node = (ASTNode_t*)StoichiometryMath_getMath(SpeciesReference_getStoichiometryMath(myRe[i]->products[j]->origin));
      if(node != NULL){/* l2v4 */
        node = ASTNode_deepCopy(node);
        myRe[i]->products[j]->value = 0;
        myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
      }else if(SpeciesReference_isSetStoichiometry(myRe[i]->products[j]->origin) && !SpeciesReference_isSetId(myRe[i]->products[j]->origin)){
        node = ASTNode_createWithType(AST_REAL);
        ASTNode_setReal(node, SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin));
        myRe[i]->products[j]->value = SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin);
        myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
      }else if(SpeciesReference_isSetStoichiometry(myRe[i]->products[j]->origin) && SpeciesReference_isSetId(myRe[i]->products[j]->origin)){
        node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(node, SpeciesReference_getId(myRe[i]->products[j]->origin));
        myRe[i]->products[j]->value = SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin);
        myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
      }
      if(node == NULL){
        if(my_isnan(SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin)) && SpeciesReference_isSetId(myRe[i]->products[j]->origin)){/* l3v1 */
          node = ASTNode_createWithType(AST_NAME);
          ASTNode_setName(node, SpeciesReference_getId(myRe[i]->products[j]->origin));
          myRe[i]->products[j]->value = 1.0;
          myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
        }else if(my_isnan(SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin))){/* l3v1 */
          node = ASTNode_createWithType(AST_REAL);
          ASTNode_setReal(node, 1.0);
          myRe[i]->products[j]->value = 1.0;
          myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
        }else{/* l2v4 */
          node = ASTNode_createWithType(AST_REAL);
          ASTNode_setReal(node, SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin));
          myRe[i]->products[j]->value = SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin);
          myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
        }
      }
      TRACE(("original math of %s : ", SpeciesReference_getSpecies(myRe[i]->products[j]->origin)));
      check_AST(node, NULL);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      /* unit */
      if(myRe[i]->products[j]->mySp->is_concentration){
        assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myRe[i]->products[j]->mySp->origin))), 1);
      }
      /* unit */
      /* conversion factor */
      if(Species_isSetConversionFactor(myRe[i]->products[j]->mySp->origin)){
        times_node = ASTNode_createWithType(AST_TIMES);
        conv_factor_node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(conv_factor_node, Species_getConversionFactor(myRe[i]->products[j]->mySp->origin));
        ASTNode_addChild(times_node, conv_factor_node);
        ASTNode_addChild(times_node, node);
        node = times_node;
      }else if(Model_isSetConversionFactor(m)){
        times_node = ASTNode_createWithType(AST_TIMES);
        conv_factor_node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(conv_factor_node, Model_getConversionFactor(m));
        ASTNode_addChild(times_node, conv_factor_node);
        ASTNode_addChild(times_node, node);
        node = times_node;
      }
      /* conversion factor */
      myRe[i]->products[j]->eq->math_length = get_equationf(m, myRe[i]->products[j]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem, print_interval);
      TRACE(("alterated math of %s : ", SpeciesReference_getSpecies(myRe[i]->products[j]->origin)));
      check_AST(node, NULL);
      /* ASTNode_free(node); */
      TRACE(("math of %s\n", SpeciesReference_getSpecies(myRe[i]->products[j]->origin)));
      check_math(myRe[i]->products[j]->eq);
      myRe[i]->products[j]->k[0] = 0;
      myRe[i]->products[j]->k[1] = 0;
      myRe[i]->products[j]->k[2] = 0;
      myRe[i]->products[j]->k[3] = 0;

	  /*6 stage (XXX new code)*/
	  myRe[i]->products[j]->k[4] = 0;
	  myRe[i]->products[j]->k[5] = 0;

      myRe[i]->products[j]->prev_val[0] = myRe[i]->products[j]->value;
      myRe[i]->products[j]->prev_val[1] = myRe[i]->products[j]->value;
      myRe[i]->products[j]->prev_val[2] = myRe[i]->products[j]->value;
      myRe[i]->products[j]->prev_k[0] = 0;
      myRe[i]->products[j]->prev_k[1] = 0;
      myRe[i]->products[j]->prev_k[2] = 0;
    }
    /* products fin */
    /* reactants start */
    for(j=0; j<myRe[i]->num_of_reactants; j++){
      myRe[i]->reactants[j]->depending_rule = NULL;
      node = (ASTNode_t*)StoichiometryMath_getMath(SpeciesReference_getStoichiometryMath(myRe[i]->reactants[j]->origin));
      if(node != NULL){/* l2v4 */
        node = ASTNode_deepCopy(node);
        myRe[i]->reactants[j]->value = 0;
        myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
      }else if(SpeciesReference_isSetStoichiometry(myRe[i]->reactants[j]->origin) && !SpeciesReference_isSetId(myRe[i]->reactants[j]->origin)){
        node = ASTNode_createWithType(AST_REAL);
        ASTNode_setReal(node, SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin));
        myRe[i]->reactants[j]->value = SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin);
        myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
      }else if(SpeciesReference_isSetStoichiometry(myRe[i]->reactants[j]->origin) && SpeciesReference_isSetId(myRe[i]->reactants[j]->origin)){
        node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(node, SpeciesReference_getId(myRe[i]->reactants[j]->origin));
        myRe[i]->reactants[j]->value = SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin);
        myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
      }
      if(node == NULL){
        if(my_isnan(SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin)) && SpeciesReference_isSetId(myRe[i]->reactants[j]->origin)){/* l3v1 */
          node = ASTNode_createWithType(AST_NAME);
          ASTNode_setName(node, SpeciesReference_getId(myRe[i]->reactants[j]->origin));
          myRe[i]->reactants[j]->value = 1.0;
          myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
        }else if(my_isnan(SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin))){/* l3v1 */
          node = ASTNode_createWithType(AST_REAL);
          ASTNode_setReal(node, 1.0);
          myRe[i]->reactants[j]->value = 1.0;
          myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
        }else{/* l2v4 */
          node = ASTNode_createWithType(AST_REAL);
          ASTNode_setReal(node, SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin));
          myRe[i]->reactants[j]->value = SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin);
          myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
        }
      }
      TRACE(("original math of %s : ", SpeciesReference_getSpecies(myRe[i]->reactants[j]->origin)));
      check_AST(node, NULL);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      /* unit */
      if(myRe[i]->reactants[j]->mySp->is_concentration){
        assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myRe[i]->reactants[j]->mySp->origin))), 1);
      }
      /* unit */
      /* conversion factor */
      if(Species_isSetConversionFactor(myRe[i]->reactants[j]->mySp->origin)){
        times_node = ASTNode_createWithType(AST_TIMES);
        conv_factor_node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(conv_factor_node, Species_getConversionFactor(myRe[i]->reactants[j]->mySp->origin));
        ASTNode_addChild(times_node, conv_factor_node);
        ASTNode_addChild(times_node, node);
        node = times_node;
      }else if(Model_isSetConversionFactor(m)){
        times_node = ASTNode_createWithType(AST_TIMES);
        conv_factor_node = ASTNode_createWithType(AST_NAME);
        ASTNode_setName(conv_factor_node, Model_getConversionFactor(m));
        ASTNode_addChild(times_node, conv_factor_node);
        ASTNode_addChild(times_node, node);
        node = times_node;
      }
      /* conversion factor */
      myRe[i]->reactants[j]->eq->math_length = get_equationf(m, myRe[i]->reactants[j]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem, print_interval);
      TRACE(("altered math of %s : ", SpeciesReference_getSpecies(myRe[i]->reactants[j]->origin)));
      check_AST(node, NULL);
      /* ASTNode_free(node); */
      TRACE(("math of %s\n", SpeciesReference_getSpecies(myRe[i]->reactants[j]->origin)));
      check_math(myRe[i]->reactants[j]->eq);
      myRe[i]->reactants[j]->k[0] = 0;      
      myRe[i]->reactants[j]->k[1] = 0;      
      myRe[i]->reactants[j]->k[2] = 0;      
      myRe[i]->reactants[j]->k[3] = 0;      

	  /*6 stage (XXX new code)*/
	  myRe[i]->reactants[j]->k[4] = 0;
	  myRe[i]->reactants[j]->k[5] = 0;

      myRe[i]->reactants[j]->prev_val[0] = myRe[i]->reactants[j]->value;
      myRe[i]->reactants[j]->prev_val[1] = myRe[i]->reactants[j]->value;
      myRe[i]->reactants[j]->prev_val[2] = myRe[i]->reactants[j]->value;
      myRe[i]->reactants[j]->prev_k[0] = 0;
      myRe[i]->reactants[j]->prev_k[1] = 0;
      myRe[i]->reactants[j]->prev_k[2] = 0;
      /* reactants fin */
    }
  }

  /* prepare reversible fast reaction */
  prepare_reversible_fast_reaction(m, myRe, mySp, myParam, myComp, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem, cp_AST);

  /* create myRules */
  for (i = 0; i < num_of_rules; i++) {
    myRu[i] = myRule_create();
    myRule_initWithModel(myRu[i], m, i);
    rule = myRule_getOrigin(myRu[i]);
    if (Rule_isRate(rule) || Rule_isAssignment(rule)) {
      myRule_initTarget(myRu[i], mySp, num_of_species, myParam, num_of_parameters, myComp, num_of_compartments, myRe, num_of_reactions);
      node = (ASTNode_t*)Rule_getMath(rule);
      node = ASTNode_deepCopy(node);
      TRACE(("original math : "));
      check_AST(node, NULL);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      /* unit */
      if(myRu[i]->target_species != NULL){
        if(myRu[i]->target_species->is_amount
            && !myRu[i]->target_species->has_only_substance_units
            && Compartment_getSpatialDimensions(myRu[i]->target_species->locating_compartment->origin) != 0){
          assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myRu[i]->target_species->origin))), 0);
        }else if(myRu[i]->target_species->is_concentration
            && myRu[i]->target_species->has_only_substance_units){
          assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myRu[i]->target_species->origin))), 1);
        }
      }
      /* unit */
      TRACE(("altered math : "));
      check_AST(node, NULL);
      myRu[i]->eq->math_length = get_equationf(m, myRu[i]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem, print_interval);
      /* ASTNode_free(node); */
      TRACE(("math\n"));
      check_math(myRule_getEquation(myRu[i]));
    }
  }

  /* create myEvents & myEventAssignments */
  for(i=0; i<num_of_events; i++){
    myEv[i] = myEvent_create();
    myEvent_initWithModel(myEv[i], m, i);
    event = myEvent_getOrigin(myEv[i]);
    node = (ASTNode_t*)Trigger_getMath(Event_getTrigger(event));
    node = ASTNode_deepCopy(node);
    TRACE(("original math of %s : ", Event_getId(myEv[i]->origin)));
    check_AST(node, NULL);
    /* alter_tree_structure(m, &node, cp_AST); */
    alter_tree_structure(m, &node, NULL, 0, cp_AST);
    myEv[i]->eq->math_length = get_equationf(m, myEv[i]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem, print_interval);
    TRACE(("altered math of %s : ", Event_getId(myEv[i]->origin)));
    check_AST(node, NULL);
    /* ASTNode_free(node); */
    TRACE(("math of %s\n", Event_getId(myEv[i]->origin)));
    check_math(myEv[i]->eq);
    myEv[i]->assignments = (myEventAssignment**)malloc(sizeof(myEventAssignment*)*Event_getNumEventAssignments(event));
    for(j=0; j<Event_getNumEventAssignments(event); j++){
      myEvAssign = myEventAssignment_create();
      myEventAssignment_initWithEvent(myEvAssign, event, j);
      myEventAssignment_initTarget(myEvAssign, mySp, num_of_species, myParam, num_of_parameters, myComp, num_of_compartments, myRe, num_of_reactions);
      myEv[i]->assignments[j] = myEvAssign;
      node = (ASTNode_t*)EventAssignment_getMath(myEv[i]->assignments[j]->origin);
      node = ASTNode_deepCopy(node);
      TRACE(("original math : "));
      check_AST(node, NULL);
      if(Event_getDelay(myEv[i]->origin) != NULL && Event_getUseValuesFromTriggerTime(myEv[i]->origin)){
        pre_ev_alter_tree_structure(&node, NULL, 0, (ASTNode_t*)Delay_getMath(Event_getDelay(myEv[i]->origin)));
        TRACE(("after pre ev alter math : "));
        check_AST(node, NULL);
        ev_alter_tree_structure(m, &node, NULL, 0, cp_AST);
        TRACE(("after ev alter math : "));
        check_AST(node, NULL);
        post_ev_alter_tree_structure(m, &node, NULL, 0);
        TRACE(("after post ev alter math : "));
        check_AST(node, NULL);
      }else{
        alter_tree_structure(m, &node, NULL, 0, cp_AST);
      }
      /* unit */
      if(myEv[i]->assignments[j]->target_species != NULL){
        if(myEv[i]->assignments[j]->target_species->is_amount
            && !myEv[i]->assignments[j]->target_species->has_only_substance_units
            && Compartment_getSpatialDimensions(myEv[i]->assignments[j]->target_species->locating_compartment->origin) != 0){
          assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myEv[i]->assignments[j]->target_species->origin))), 0);
        }else if(myEv[i]->assignments[j]->target_species->is_concentration
            && myEv[i]->assignments[j]->target_species->has_only_substance_units){
          assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myEv[i]->assignments[j]->target_species->origin))), 1);
        }
      }
      /* unit */
      TRACE(("altered math : "));
      check_AST(node, NULL);
      myEv[i]->assignments[j]->eq->math_length = get_equationf(m, myEv[i]->assignments[j]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem, print_interval);
      /* ASTNode_free(node); */
      TRACE(("math\n"));
      check_math(myEv[i]->assignments[j]->eq);
    }
    myEv[i]->event_delay = NULL;
    if(Event_getDelay(event) != NULL){
      myEv[i]->event_delay = (myDelay*)malloc(sizeof(myDelay));
      myEv[i]->event_delay->origin = Event_getDelay(event);
      node = (ASTNode_t*)Delay_getMath(myEv[i]->event_delay->origin);
      node = ASTNode_deepCopy(node);
      TRACE(("original math : "));
      check_AST(node, NULL);
      /* alter_tree_structure(m, &node, cp_AST); */
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      myEv[i]->event_delay->eq = equation_create();
      myEv[i]->event_delay->eq->math_length = get_equationf(m, myEv[i]->event_delay->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem, print_interval);
      TRACE(("altered math : "));
      check_AST(node, NULL);
      /* ASTNode_free(node); */
      TRACE(("math\n"));
      check_math(myEv[i]->event_delay->eq);
      myEv[i]->firing_times = (double*)malloc(sizeof(double)*(int)(sim_time/dt));
      for(j=0; j<(unsigned int)(sim_time/dt); j++){
        myEv[i]->firing_times[j] = sim_time+1;
      }
      myEv[i]->num_of_delayed_events_que = 0;
      myEv[i]->next_firing_index = 0;
    }
    myEv[i]->priority_eq = NULL;
    if(Event_isSetPriority(myEv[i]->origin)){
      myEv[i]->priority_eq = equation_create();
      node = (ASTNode_t*)Priority_getMath(Event_getPriority(myEv[i]->origin));
      node = ASTNode_deepCopy(node);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      myEv[i]->priority_eq->math_length = get_equationf(m, myEv[i]->priority_eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem, print_interval);
    }
  }

  /* prepare Algebraic */
  num_of_algebraic_rules = 0;
  for(i=0; i<num_of_rules; i++){
    if(myRu[i]->is_algebraic){
      num_of_algebraic_rules++;
    }
  }
  if(num_of_algebraic_rules != 0){
    *myAlgEq = (myAlgebraicEquations*)malloc(sizeof(myAlgebraicEquations));
    algEq = *myAlgEq;
    algEq->num_of_algebraic_rules = num_of_algebraic_rules;
    algEq->coefficient_matrix = NULL;
    algEq->constant_vector = NULL;
    algEq->coefficient = NULL;
    algEq->constant = NULL;
    algEq->target_species = NULL;
    algEq->target_parameter = NULL;
    algEq->target_compartment = NULL;
    for(i=0; i<MAX_ALGEBRAIC_VARIABLES; i++){
      algEq->alg_target_species[i] = NULL;
      algEq->alg_target_parameter[i] = NULL;
      algEq->alg_target_compartment[i] = NULL;
    }
    algEq->num_of_alg_target_sp = 0;
    algEq->num_of_alg_target_param = 0;
    algEq->num_of_alg_target_comp = 0;
    algEq->num_of_algebraic_variables = 0;
    if(algEq->num_of_algebraic_rules > 1){
      algEq->coefficient_matrix = (equation***)malloc(sizeof(equation**)*algEq->num_of_algebraic_rules);
      for(i=0; i<algEq->num_of_algebraic_rules; i++){
        algEq->coefficient_matrix[i] = (equation**)malloc(sizeof(equation*)*algEq->num_of_algebraic_rules);
        for(j=0; j<algEq->num_of_algebraic_rules; j++){
          algEq->coefficient_matrix[i][j] = equation_create();
          algEq->coefficient_matrix[i][j]->math_length = 0;
        }
      }
      algEq->constant_vector = (equation**)malloc(sizeof(equation*)*algEq->num_of_algebraic_rules);
      for(i=0; i<algEq->num_of_algebraic_rules; i++){
        algEq->constant_vector[i] = equation_create();
      }
    }else{
      algEq->coefficient = equation_create();
      algEq->constant = equation_create();
    }
    TRACE(("prepare algebraic start\n"));
    prepare_algebraic(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, algEq, sim_time, dt, time, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem, cp_AST);
    TRACE(("prepare algebraic finish\n"));
    if(algEq->num_of_algebraic_rules > 1){
      for(i=0; i<algEq->num_of_algebraic_rules; i++){
        for(j=0; j<algEq->num_of_algebraic_rules; j++){
          TRACE(("math of coefficient matrix[%d][%d] is\n", i, j));
          check_math(algEq->coefficient_matrix[i][j]);
        }
      }
      for(i=0; i<algEq->num_of_algebraic_rules; i++){
        TRACE(("constant vector[%d] is\n", i));
        check_math(algEq->constant_vector[i]);
      }
    }else{
      TRACE(("math of coefficient is\n"));
      check_math(algEq->coefficient);
      TRACE(("math of constant is\n"));
      check_math(algEq->constant);
    }
    if(algEq->num_of_algebraic_rules > 1){
      flag = 0;
      for(i=0; i<num_of_species; i++){
        for(j=0; j<algEq->num_of_algebraic_variables; j++){
          if(strcmp(algEq->variables_id[j], Species_getId(mySp[i]->origin)) == 0){
            algEq->alg_target_species[algEq->num_of_alg_target_sp] = (myAlgTargetSp*)malloc(sizeof(myAlgTargetSp));
            algEq->alg_target_species[algEq->num_of_alg_target_sp]->target_species = mySp[i];
            algEq->alg_target_species[algEq->num_of_alg_target_sp]->order = j;
            algEq->num_of_alg_target_sp++;
            flag = 1;
            break;
          }
        }
      }
      if(!flag){
        for(i=0; i<num_of_parameters; i++){
          for(j=0; j<algEq->num_of_algebraic_variables; j++){
            if(strcmp(algEq->variables_id[j], Parameter_getId(myParam[i]->origin)) == 0){
              algEq->alg_target_parameter[algEq->num_of_alg_target_param] = (myAlgTargetParam*)malloc(sizeof(myAlgTargetParam));
              algEq->alg_target_parameter[algEq->num_of_alg_target_param]->target_parameter = myParam[i];
              algEq->alg_target_parameter[algEq->num_of_alg_target_param]->order = j;
              algEq->num_of_alg_target_param++;
              flag = 1;
              break;
            }
          }
        }
      }
      if(!flag){
        for(i=0; i<num_of_compartments; i++){
          for(j=0; j<algEq->num_of_algebraic_variables; j++){
            if(strcmp(algEq->variables_id[j], Compartment_getId(myComp[i]->origin)) == 0){
              algEq->alg_target_compartment[algEq->num_of_alg_target_comp] = (myAlgTargetComp*)malloc(sizeof(myAlgTargetComp));
              algEq->alg_target_compartment[algEq->num_of_alg_target_comp]->target_compartment = myComp[i];
              algEq->alg_target_compartment[algEq->num_of_alg_target_comp]->order = j;
              algEq->num_of_alg_target_comp++;
              flag = 1;
              break;
            }
          }
        }
      }
    }else{
      flag = 0;
      for(i=0; i<num_of_species; i++){
        if(strcmp(algEq->variables_id[0], Species_getId(mySp[i]->origin)) == 0){
          algEq->target_species = mySp[i];
          flag = 1;
          break;
        }
      }
      if(!flag){
        for(i=0; i<num_of_parameters; i++){
          if(strcmp(algEq->variables_id[0], Parameter_getId(myParam[i]->origin)) == 0){
            algEq->target_parameter = myParam[i];
            flag = 1;
            break;
          }
        }
      }
      if(!flag){
        for(i=0; i<num_of_compartments; i++){
          if(strcmp(algEq->variables_id[0], Compartment_getId(myComp[i]->origin)) == 0){
            algEq->target_compartment = myComp[i];
            flag = 1;
            break;
          }
        }
      }
    }
  }
}

void free_mySBML_objects(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], myReaction *myRe[], myRule *myRu[], myEvent *myEv[], myInitialAssignment *myInitAssign[], myAlgebraicEquations *myAlgEq, timeVariantAssignments *timeVarAssign, allocated_memory *mem, copied_AST *cp_AST){
  unsigned int i, j;

  for (i = 0; i < Model_getNumSpecies(m); i++) {
    mySpecies_free(mySp[i]);
  }
  free(mySp);

  for (i = 0; i < Model_getNumParameters(m); i++) {
    myParameter_free(myParam[i]);
  }
  free(myParam);

  for (i = 0; i < Model_getNumCompartments(m); i++) {
    myCompartment_free(myComp[i]);
  }
  free(myComp);

  for (i = 0; i < Model_getNumReactions(m); i++) {
    myReaction_free(myRe[i]);
  }
  free(myRe);

  for (i = 0; i < Model_getNumRules(m); i++) {
    myRule_free(myRu[i]);
  }
  free(myRu);

  for (i = 0; i < Model_getNumEvents(m); i++) {
    myEvent_free(myEv[i]);
  }
  free(myEv);

  for (i = 0; i < Model_getNumInitialAssignments(m); i++) {
    myInitialAssignment_free(myInitAssign[i]);
  }
  free(myInitAssign);

  if(myAlgEq != NULL){
    if(myAlgEq->num_of_algebraic_variables > 1){
      for(i=0; i<myAlgEq->num_of_algebraic_variables; i++){
        for(j=0; j<myAlgEq->num_of_algebraic_variables; j++){
          free(myAlgEq->coefficient_matrix[i][j]);
        }
        free(myAlgEq->coefficient_matrix[i]);
      }
      free(myAlgEq->coefficient_matrix);
      for(i=0; i<myAlgEq->num_of_algebraic_variables; i++){
        free(myAlgEq->constant_vector[i]);
      }
      free(myAlgEq->constant_vector);
      for(i=0; i<myAlgEq->num_of_alg_target_sp; i++){
        free(myAlgEq->alg_target_species[i]);
      }
      for(i=0; i<myAlgEq->num_of_alg_target_param; i++){
        free(myAlgEq->alg_target_parameter[i]);
      }
      for(i=0; i<myAlgEq->num_of_alg_target_comp; i++){
        free(myAlgEq->alg_target_compartment[i]);
      }
    }else{
      free(myAlgEq->coefficient);
      free(myAlgEq->constant);
    }
    free(myAlgEq);
  }
  for(i=0; i<timeVarAssign->num_of_time_variant_assignments; i++){
    free(timeVarAssign->eq[i]);
  }
  free(timeVarAssign);

  allocated_memory_free(mem);
  copied_AST_free(cp_AST);

  TRACE(("all allocated memory is freeed\n"));
}

void realloc_mySBML_objects(Model_t *m, mySpecies *sp[], unsigned int num_of_species, myParameter *param[], unsigned int num_of_parameters, myCompartment *comp[], unsigned int num_of_compartments, myReaction *re[], unsigned int num_of_reactions, myRule *rule[], myEvent *ev[], unsigned int num_of_events, myInitialAssignment *myInitAssign[], timeVariantAssignments **timeVarAssign, copied_AST *cp_AST, double sim_time, int max_index) {
	unsigned int i, j, l;
	int k;
	int new_max_index = (max_index - 1) * 2 + 1;
	/* for substitution of reallocated address */
	int flag = 0;
	ASTNode_t* node;
	ASTNode_t *times_node, *conv_factor_node;
	unsigned int math_length = 0;
	/* 1) reallocate delay value arrays */
	/* species */
	for (i = 0; i < num_of_species; i++) {
		if (sp[i]->delay_val != NULL) {
      mySpecies_reallocDelayVal(sp[i], new_max_index, 6);
			flag = 1;
		}
	}
	/* parameter */
	for(i=0; i<num_of_parameters; i++){
		if (param[i]->delay_val != NULL) {
      myParameter_reallocDelayVal(param[i], new_max_index, 6);
			flag = 1;
		}
	}
	/* compartment*/
	for (i = 0; i < num_of_compartments; i++) {
		if (comp[i]->delay_val != NULL) {
      myCompartment_reallocDelayVal(comp[i], new_max_index, 6);
			flag = 1;
		}
	}
	/* products and reactants */
	for(i=0; i<num_of_reactions; i++){
		for(j=0; j<re[i]->num_of_products; j++){
			if(re[i]->products[j]->delay_val != NULL){
				re[i]->products[j]->delay_val = (double**)realloc(re[i]->products[j]->delay_val, sizeof(double*)* new_max_index);
				for(k=0; k<=new_max_index; k++) {
					if (k <max_index) {
						re[i]->products[j]->delay_val[k] = (double *)realloc(re[i]->products[j]->delay_val[k], sizeof(double) * 6);
					}else {
						re[i]->products[j]->delay_val[k] = (double *)malloc(sizeof(double) * 6);
					}
				}
				flag = 1;
			}
		}
		for(j=0; j<re[i]->num_of_reactants; j++){
			if(re[i]->reactants[j]->delay_val != NULL){
				re[i]->reactants[j]->delay_val = (double**)realloc(re[i]->reactants[j]->delay_val, sizeof(double*)* new_max_index);
				for(k=0; k<=new_max_index; k++) {
					if (k <max_index) {
						re[i]->reactants[j]->delay_val[k] = (double *)realloc(re[i]->reactants[j]->delay_val[k], sizeof(double) * 6);
					}else {
						re[i]->reactants[j]->delay_val[k] = (double *)malloc(sizeof(double) * 6);
					}
				}
				flag = 1;
			}
		}
	}
	if(flag) {
		/* connect delay_val with all the equation->delay_number */
		/* initial assignments */
		for(i=0; i<Model_getNumInitialAssignments(m); i++){
			node = (ASTNode_t*)InitialAssignment_getMath(myInitAssign[i]->origin);
			node = ASTNode_deepCopy(node);
			check_AST(node, NULL);
			/* alter_tree_structure */
			alter_tree_structure(m, &node, NULL, 0, cp_AST);
			math_length = connect_delayval_with_eq(m, myInitAssign[i]->eq, sp, param, comp, re, node, 0);
		}
		/* time var assignments */
		(*timeVarAssign)->num_of_time_variant_assignments = 0;
		for(i=0; i<Model_getNumRules(m); i++){
			node = (ASTNode_t*)Rule_getMath(rule[i]->origin);
			if(Rule_isAssignment(rule[i]->origin) && include_time(node, 0)){
				node = ASTNode_deepCopy(node);
				check_AST(node, NULL);
				alter_tree_structure(m, &node, NULL, 0, cp_AST);
				for(j=0; j<Model_getNumSpecies(m); j++){
					if(strcmp(Species_getId(sp[j]->origin), (char*)Rule_getVariable(rule[i]->origin)) == 0){
						if(sp[j]->is_amount
						   && !sp[j]->has_only_substance_units
						   && Compartment_getSpatialDimensions(sp[j]->locating_compartment->origin) != 0){
							assignment_alter_tree_structure(&node, (char*)Compartment_getId(sp[j]->locating_compartment->origin), 0);
						}else if(sp[j]->is_concentration
								 && sp[j]->has_only_substance_units){
							assignment_alter_tree_structure(&node, (char*)Compartment_getId(sp[j]->locating_compartment->origin), 1);
						}
						break;
					}
				}
				check_AST(node, NULL);
				math_length = connect_delayval_with_eq(m, (*timeVarAssign)->eq[(*timeVarAssign)->num_of_time_variant_assignments], sp, param, comp, re, node, 0);
			}
		}
		/* reaction equation */
		for(i=0; i<num_of_reactions; i++){
			node = (ASTNode_t*)KineticLaw_getMath(Reaction_getKineticLaw(re[i]->origin));
			node = ASTNode_deepCopy(node);
			check_AST(node, NULL);
			alter_tree_structure(m, &node, NULL, 0, cp_AST);
			set_local_para_as_value(node, Reaction_getKineticLaw(re[i]->origin));
			check_AST(node, NULL);
			math_length = connect_delayval_with_eq(m, re[i]->eq, sp, param, comp, re, node, 0);
			/* products */
			for(j=0; j<re[i]->num_of_products; j++){
				node = (ASTNode_t*)StoichiometryMath_getMath(SpeciesReference_getStoichiometryMath(re[i]->products[j]->origin));
				if(node != NULL){/* l2v4 */
					node = ASTNode_deepCopy(node);
				}else if(SpeciesReference_isSetStoichiometry(re[i]->products[j]->origin) && !SpeciesReference_isSetId(re[i]->products[j]->origin)){
					node = ASTNode_createWithType(AST_REAL);
					ASTNode_setReal(node, SpeciesReference_getStoichiometry(re[i]->products[j]->origin));
				}else if(SpeciesReference_isSetStoichiometry(re[i]->products[j]->origin) && SpeciesReference_isSetId(re[i]->products[j]->origin)){
					node = ASTNode_createWithType(AST_NAME);
					ASTNode_setName(node, SpeciesReference_getId(re[i]->products[j]->origin));
				}
				if(node == NULL){
					if(my_isnan(SpeciesReference_getStoichiometry(re[i]->products[j]->origin)) && SpeciesReference_isSetId(re[i]->products[j]->origin)){/* l3v1 */
						node = ASTNode_createWithType(AST_NAME);
						ASTNode_setName(node, SpeciesReference_getId(re[i]->products[j]->origin));
					}else if(my_isnan(SpeciesReference_getStoichiometry(re[i]->products[j]->origin))){/* l3v1 */
						node = ASTNode_createWithType(AST_REAL);
						ASTNode_setReal(node, 1.0);
					}else{/* l2v4 */
						node = ASTNode_createWithType(AST_REAL);
						ASTNode_setReal(node, SpeciesReference_getStoichiometry(re[i]->products[j]->origin));
					}
				}
				check_AST(node, NULL);
				alter_tree_structure(m, &node, NULL, 0, cp_AST);
				/* unit */
				if(re[i]->products[j]->mySp->is_concentration){
					assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(re[i]->products[j]->mySp->origin))), 1);
				}
				/* unit */
				/* conversion factor */
				if(Species_isSetConversionFactor(re[i]->products[j]->mySp->origin)){
					times_node = ASTNode_createWithType(AST_TIMES);
					conv_factor_node = ASTNode_createWithType(AST_NAME);
					ASTNode_setName(conv_factor_node, Species_getConversionFactor(re[i]->products[j]->mySp->origin));
					ASTNode_addChild(times_node, conv_factor_node);
					ASTNode_addChild(times_node, node);
					node = times_node;
				}else if(Model_isSetConversionFactor(m)){
					times_node = ASTNode_createWithType(AST_TIMES);
					conv_factor_node = ASTNode_createWithType(AST_NAME);
					ASTNode_setName(conv_factor_node, Model_getConversionFactor(m));
					ASTNode_addChild(times_node, conv_factor_node);
					ASTNode_addChild(times_node, node);
					node = times_node;
				}
				/* conversion factor */
				math_length = connect_delayval_with_eq(m, re[i]->products[j]->eq, sp, param, comp, re, node, 0);
				check_AST(node, NULL);
			}
			/* reactants */
			for(j=0; j<re[i]->num_of_reactants; j++){
				node = (ASTNode_t*)StoichiometryMath_getMath(SpeciesReference_getStoichiometryMath(re[i]->reactants[j]->origin));
				if(node != NULL){/* l2v4 */
					node = ASTNode_deepCopy(node);
				}else if(SpeciesReference_isSetStoichiometry(re[i]->reactants[j]->origin) && !SpeciesReference_isSetId(re[i]->reactants[j]->origin)){
					node = ASTNode_createWithType(AST_REAL);
					ASTNode_setReal(node, SpeciesReference_getStoichiometry(re[i]->reactants[j]->origin));
				}else if(SpeciesReference_isSetStoichiometry(re[i]->reactants[j]->origin) && SpeciesReference_isSetId(re[i]->reactants[j]->origin)){
					node = ASTNode_createWithType(AST_NAME);
					ASTNode_setName(node, SpeciesReference_getId(re[i]->reactants[j]->origin));
				}
				if(node == NULL){
					if(my_isnan(SpeciesReference_getStoichiometry(re[i]->reactants[j]->origin)) && SpeciesReference_isSetId(re[i]->reactants[j]->origin)){/* l3v1 */
						node = ASTNode_createWithType(AST_NAME);
						ASTNode_setName(node, SpeciesReference_getId(re[i]->reactants[j]->origin));
					}else if(my_isnan(SpeciesReference_getStoichiometry(re[i]->reactants[j]->origin))){/* l3v1 */
						node = ASTNode_createWithType(AST_REAL);
						ASTNode_setReal(node, 1.0);
					}else{/* l2v4 */
						node = ASTNode_createWithType(AST_REAL);
						ASTNode_setReal(node, SpeciesReference_getStoichiometry(re[i]->reactants[j]->origin));
					}
				}
				check_AST(node, NULL);
				alter_tree_structure(m, &node, NULL, 0, cp_AST);
				/* unit */
				if(re[i]->reactants[j]->mySp->is_concentration){
					assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(re[i]->reactants[j]->mySp->origin))), 1);
				}
				/* unit */
				/* conversion factor */
				if(Species_isSetConversionFactor(re[i]->reactants[j]->mySp->origin)){
					times_node = ASTNode_createWithType(AST_TIMES);
					conv_factor_node = ASTNode_createWithType(AST_NAME);
					ASTNode_setName(conv_factor_node, Species_getConversionFactor(re[i]->reactants[j]->mySp->origin));
					ASTNode_addChild(times_node, conv_factor_node);
					ASTNode_addChild(times_node, node);
					node = times_node;
				}else if(Model_isSetConversionFactor(m)){
					times_node = ASTNode_createWithType(AST_TIMES);
					conv_factor_node = ASTNode_createWithType(AST_NAME);
					ASTNode_setName(conv_factor_node, Model_getConversionFactor(m));
					ASTNode_addChild(times_node, conv_factor_node);
					ASTNode_addChild(times_node, node);
					node = times_node;
				}
				/* conversion factor */
				math_length =connect_delayval_with_eq(m, re[i]->reactants[j]->eq, sp, param, comp, re, node, 0);
				check_AST(node, NULL);
			}
		}

		/* create myRules */
		for(i=0; i<Model_getNumRules(m); i++){
			if(Rule_isRate(rule[i]->origin) || Rule_isAssignment(rule[i]->origin)){
				flag = 1;
				for(j=0; j<num_of_species; j++){
					if(strcmp(Rule_getVariable(rule[i]->origin), Species_getId(sp[j]->origin)) == 0){
						flag = 0;
						break;
					}
				}
				if(flag){
					for(j=0; j<num_of_parameters; j++){
						if(strcmp(Rule_getVariable(rule[i]->origin), Parameter_getId(param[j]->origin)) == 0){
							flag = 0;
							break;
						}
					}
				}
				if(flag){
					for(j=0; j<num_of_compartments; j++){
						if(strcmp(Rule_getVariable(rule[i]->origin), Compartment_getId(comp[j]->origin)) == 0){
							flag = 0;
							break;
						}
					}
				}
				if(flag){
					for(j=0; j<num_of_reactions; j++){
						for(l=0; l<re[j]->num_of_products; l++){
							if(SpeciesReference_isSetId(re[j]->products[l]->origin)
							   && strcmp(Rule_getVariable(rule[i]->origin), SpeciesReference_getId(re[j]->products[l]->origin)) == 0){
								flag = 0;
								break;
							}
						}
						if(!flag){
							break;
						}
						for(l=0; l<re[j]->num_of_reactants; l++){
							if(SpeciesReference_isSetId(re[j]->reactants[l]->origin)
							   && strcmp(Rule_getVariable(rule[i]->origin), SpeciesReference_getId(re[j]->reactants[l]->origin)) == 0){
								flag = 0;
								break;
							}
						}
						if(!flag){
							break;
						}
					}
				}
				node = (ASTNode_t*)Rule_getMath(rule[i]->origin);
				node = ASTNode_deepCopy(node);
				TRACE(("original math : "));
				check_AST(node, NULL);
				alter_tree_structure(m, &node, NULL, 0, cp_AST);
				/* unit */
				if(rule[i]->target_species != NULL){
					if(rule[i]->target_species->is_amount
					   && !rule[i]->target_species->has_only_substance_units
					   && Compartment_getSpatialDimensions(rule[i]->target_species->locating_compartment->origin) != 0){
						assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(rule[i]->target_species->origin))), 0);
					}else if(rule[i]->target_species->is_concentration
							 && rule[i]->target_species->has_only_substance_units){
						assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(rule[i]->target_species->origin))), 1);
					}
				}
				/* unit */
				TRACE(("altered math : "));
				check_AST(node, NULL);
				math_length = connect_delayval_with_eq(m, rule[i]->eq, sp, param, comp, re, node, 0);
				/* ASTNode_free(node); */
				TRACE(("math\n"));
				check_math(rule[i]->eq);
			}
		}

		/* event */
		for(i=0; i<num_of_events; i++){
			node = (ASTNode_t*)Trigger_getMath(Event_getTrigger(ev[i]->origin));
			node = ASTNode_deepCopy(node);
			check_AST(node, NULL);
			alter_tree_structure(m, &node, NULL, 0, cp_AST);
			math_length = connect_delayval_with_eq(m, ev[i]->eq, sp, param, comp, re, node, 0);
			check_AST(node, NULL);
			if (Event_getDelay(ev[i]->origin) != NULL){
				node = (ASTNode_t*)Delay_getMath(ev[i]->event_delay->origin);
				node = ASTNode_deepCopy(node);
				check_AST(node, NULL);
				alter_tree_structure(m, &node, NULL, 0, cp_AST);
				math_length = connect_delayval_with_eq(m, ev[i]->event_delay->eq, sp, param, comp, re, node, 0);
			}
			for(j=0; j<Event_getNumEventAssignments(ev[i]->origin); j++){
				node = (ASTNode_t*)EventAssignment_getMath(ev[i]->assignments[j]->origin);
				node = ASTNode_deepCopy(node);
				check_AST(node, NULL);
				if(Event_getDelay(ev[i]->origin) != NULL && Event_getUseValuesFromTriggerTime(ev[i]->origin)){
					pre_ev_alter_tree_structure(&node, NULL, 0, (ASTNode_t*)Delay_getMath(Event_getDelay(ev[i]->origin)));
					TRACE(("after pre ev alter math : "));
					check_AST(node, NULL);
					ev_alter_tree_structure(m, &node, NULL, 0, cp_AST);
					TRACE(("after ev alter math : "));
					check_AST(node, NULL);
					post_ev_alter_tree_structure(m, &node, NULL, 0);
					TRACE(("after post ev alter math : "));
					check_AST(node, NULL);
				}else{
					alter_tree_structure(m, &node, NULL, 0, cp_AST);
				}
				/* unit */
				if(ev[i]->assignments[j]->target_species != NULL){
					if(ev[i]->assignments[j]->target_species->is_amount
					   && !ev[i]->assignments[j]->target_species->has_only_substance_units
					   && Compartment_getSpatialDimensions(ev[i]->assignments[j]->target_species->locating_compartment->origin) != 0){
						assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(ev[i]->assignments[j]->target_species->origin))), 0);
					}else if(ev[i]->assignments[j]->target_species->is_concentration
							 && ev[i]->assignments[j]->target_species->has_only_substance_units){
						assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(ev[i]->assignments[j]->target_species->origin))), 1);
					}
				}
				/* unit */
				check_AST(node, NULL);
				math_length = connect_delayval_with_eq(m, ev[i]->assignments[j]->eq, sp, param, comp, re, node, 0);
			}
		}
	}
	/* 2) reallocate event firing times
	   (note that new max index of array is smaller than new_max_index by 1)*/
	for(i=0; i<num_of_events; i++) {
		if(ev[i] != NULL && Event_getDelay(ev[i]->origin) != NULL){
			ev[i]->firing_times = (double*)realloc(ev[i]->firing_times, sizeof(double)*(int)new_max_index-1);
			for(k=max_index; k<new_max_index-1; k++){
				ev[i]->firing_times[k] = sim_time+1;
			}
		}
	}
}

void reallocate_result_objects(myResult* res, double** value_time_p_fordelay, unsigned  int max_result_index){
	double* head_time_e;

	/* save the value of max_index in Result */
	res->num_of_delay_rows = max_result_index * 2;

	/* reallocate value arrays in Result */
	head_time_e = (double *)realloc(res->values_time_fordelay, sizeof(double) * res->num_of_delay_rows);

	/* error check */
	if (!head_time_e) {
		fprintf(stderr, "failed in reallocation of time memory.\n");
		exit(1);
	}else{
		res->values_time_fordelay = head_time_e;
		head_time_e += max_result_index;
		*value_time_p_fordelay = head_time_e;
	}
}

unsigned int connect_delayval_with_eq(Model_t *m, equation *eq, mySpecies *sp[], myParameter *param[], myCompartment *comp[], myReaction *re[], ASTNode_t *node, int index){
	unsigned int i,j;
	int flag;
	const char* name;
	ASTNode_t *left, *right, *comp_node;
	if((ASTNode_getType(node) == AST_LOGICAL_AND
        || ASTNode_getType(node) == AST_LOGICAL_OR
        || ASTNode_getType(node) == AST_LOGICAL_XOR)
	   && ASTNode_getNumChildren(node) > 2){
		ASTNode_reduceToBinary(node);
	}
	if(ASTNode_getType(node) == AST_FUNCTION_DELAY){
		left = ASTNode_getLeftChild(node);
		comp_node = NULL;
		if(ASTNode_getType(left) != AST_NAME){
			comp_node = ASTNode_getRightChild(left);
			left = ASTNode_getLeftChild(left);
		}
		name = ASTNode_getName(left);
		flag = 1;
		/* species */
		for(i=0; i<Model_getNumSpecies(m); i++){
			if(strcmp(name, Species_getId(sp[i]->origin)) == 0){
				/* connection */
				eq->delay_number[index] = sp[i]->delay_val;
				if(comp_node != NULL){
					for(j=0; j<Model_getNumCompartments(m); j++){
						if(strcmp(ASTNode_getName(comp_node), Compartment_getId(comp[j]->origin)) == 0){
							/* connection */
							eq->delay_comp_size[index] = comp[j]->delay_val;
							break;
						}
					}
				}
				index++;
				flag = 0;
				break;
			}
		}
		if(flag){
			/* parameter */
			for(i=0; i<Model_getNumParameters(m); i++){
				if(strcmp(name, Parameter_getId(param[i]->origin)) == 0){
					/* connection */
					eq->delay_number[index] = param[i]->delay_val;
					index++;
					flag = 0;
					break;
				}
			}
		}
		if(flag){
			/* compartment */
			for(i=0; i<Model_getNumCompartments(m); i++){
				if(strcmp(name, Compartment_getId(comp[i]->origin)) == 0){
					/* connection */
					eq->delay_number[index] = comp[i]->delay_val;
					index++;
					flag = 0;
					break;
				}
			}
		}
		if(flag){
			for(i=0; i<Model_getNumReactions(m); i++){
				/* products */
				for(j=0; j<re[i]->num_of_products; j++){
					if(SpeciesReference_isSetId(re[i]->products[j]->origin)
					   && strcmp(name, SpeciesReference_getId(re[i]->products[j]->origin)) == 0){
						/* connection */
						eq->delay_number[index] = re[i]->products[j]->delay_val;
						index++;
						flag = 0;
						break;
					}
				}
				if(!flag){
					break;
				}
				/* reactants */
				for(j=0; j<re[i]->num_of_reactants; j++){
					if(SpeciesReference_isSetId(re[i]->reactants[j]->origin)
					   && strcmp(name, SpeciesReference_getId(re[i]->reactants[j]->origin)) == 0){
						/* connection */
						eq->delay_number[index] = re[i]->reactants[j]->delay_val;
						index++;
						flag = 0;
						break;
					}
				}
				if(!flag){
					break;
				}
			}
		}
	}else if((left=ASTNode_getLeftChild(node)) != NULL){
		index = connect_delayval_with_eq(m, eq, sp, param, comp, re, left, index);
	}
	if((right=ASTNode_getRightChild(node)) != NULL){
		index = connect_delayval_with_eq(m, eq, sp, param, comp, re, right, index);
	}

  if(ASTNode_isOperator(node)
      || ASTNode_isFunction(node)
      || ASTNode_isBoolean(node)){
    index++;
  }else if(ASTNode_getType(node) == AST_NAME){
    name = ASTNode_getName(node);
    flag = 1;
    for(i=0; i<Model_getNumSpecies(m); i++){
      if(strcmp(name, Species_getId(sp[i]->origin)) == 0){
        index++;
        flag = 0;
        break;
      }
    }
    if(flag){
      for(i=0; i<Model_getNumParameters(m); i++){
        if(strcmp(name, Parameter_getId(param[i]->origin)) == 0){
          index++;
          flag = 0;
          break;
        }
      }
    }
    if(flag){
      for(i=0; i<Model_getNumCompartments(m); i++){
        if(strcmp(name, Compartment_getId(comp[i]->origin)) == 0){
          index++;
          flag = 0;
          break;
        }
      }
    }
    if(flag){
      for(i=0; i<Model_getNumReactions(m); i++){
        for(j=0; j<re[i]->num_of_products; j++){
          if(SpeciesReference_isSetId(re[i]->products[j]->origin)
              && strcmp(name, SpeciesReference_getId(re[i]->products[j]->origin)) == 0){
            index++;
            flag = 0;
            break;
          }
        }
        if(!flag){
          break;
        }
        for(j=0; j<re[i]->num_of_reactants; j++){
          if(SpeciesReference_isSetId(re[i]->reactants[j]->origin)
              && strcmp(name, SpeciesReference_getId(re[i]->reactants[j]->origin)) == 0){
            index++;
            flag = 0;
            break;
          }
        }
        if(!flag){
          break;
        }
      }
    }
    if(flag){
      if(strcmp(name, "time") == 0
          || strcmp(name, "t") == 0
          || strcmp(name, "s") == 0){
        index++;
      }
    }
  }else if(ASTNode_getType(node) == AST_NAME_TIME){
    index++;
  }else if(ASTNode_getType(node) == AST_NAME_AVOGADRO){
    index++;
  }else if(ASTNode_getType(node) == AST_CONSTANT_E
      || ASTNode_getType(node) == AST_CONSTANT_PI){
    index++;
  }else{
    index++;
  }
	return index;
}

static int include_time(ASTNode_t *node, int flag) {
  unsigned int i;
  char *name;

  for (i = 0; i < ASTNode_getNumChildren(node); i++) {
    flag = include_time(ASTNode_getChild(node, i), flag);
  }
  if (ASTNode_getType(node) == AST_NAME_TIME) {
    flag = 1;
  } else if (ASTNode_getType(node) == AST_NAME) {
    name = (char*)ASTNode_getName(node);
    if (strcmp(name, "time") == 0
        || strcmp(name, "t") == 0
        || strcmp(name, "s") == 0){
      flag = 1;
    }
  }
  return flag;
}

static void locate_species_in_compartment(
    mySpecies **species, int num_of_species,
    myCompartment **compartments, int num_of_compartments) {
  int i, j;
  const char *cid;
  for (i = 0; i < num_of_species; i++) {
    cid = Species_getCompartment(species[i]->origin);
    for (j = 0; j < num_of_compartments; j++) {
      if (strcmp(cid, Compartment_getId(compartments[j]->origin)) == 0) {
        species[i]->locating_compartment = compartments[j];
        compartments[j]->including_species[compartments[j]->num_of_including_species] = species[i];
        compartments[j]->num_of_including_species++;
      }
    }
  }
}

static void create_species(mySpecies *species[], Model_t *model) {
  unsigned int i;
  unsigned int num_of_species = Model_getNumSpecies(model);
  for (i = 0; i < num_of_species; i++) {
    species[i] = mySpecies_create();
    mySpecies_initWithModel(species[i], model, i);
  }
}

static void create_parameters(myParameter *parameters[], Model_t *model) {
  unsigned int i;
  unsigned int num_of_parameters = Model_getNumParameters(model);
  for (i = 0; i < num_of_parameters; i++) {
    parameters[i] = myParameter_create();
    myParameter_initWithModel(parameters[i], model, i);
  }
}

static void create_compartments(myCompartment *compartments[], Model_t *model) {
  unsigned int i;
  unsigned int num_of_compartments = Model_getNumCompartments(model);
  for (i = 0; i < num_of_compartments; i++) {
    compartments[i] = myCompartment_create();
    myCompartment_initWithModel(compartments[i], model, i);
  }
}

static void create_reactions(myReaction *reactions[], mySpecies *species[], Model_t *model) {
  unsigned int i;
  unsigned int num_of_reactions = Model_getNumReactions(model);
  unsigned int num_of_species = Model_getNumSpecies(model);
  for (i = 0; i < num_of_reactions; i++) {
    reactions[i] = myReaction_create();
    myReaction_initWithModel(reactions[i], model, i);
    myReaction_initProducts(reactions[i], species, num_of_species);
    myReaction_initReactants(reactions[i], species, num_of_species);
  }
}

