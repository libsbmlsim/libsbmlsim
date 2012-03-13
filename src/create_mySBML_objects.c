#include "libsbmlsim/libsbmlsim.h"

int include_time(ASTNode_t *node, int flag){
  int i;
  char *name;

  for(i=0; i<ASTNode_getNumChildren(node); i++){
    flag = include_time(ASTNode_getChild(node, i), flag);
  }
  if(ASTNode_getType(node) == AST_NAME_TIME){
    flag = 1;
  }else if(ASTNode_getType(node) == AST_NAME){
    name = (char*)ASTNode_getName(node);
    if(strcmp(name, "time") == 0
        || strcmp(name, "t") == 0
        || strcmp(name, "s") == 0){
      flag = 1;
    }
  }
  return flag;
}

//create my SBML obejects for efficient simulations
void create_mySBML_objects(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], myReaction *myRe[], myRule *myRu[], myEvent *myEv[], myInitialAssignment *myInitAssign[], myAlgebraicEquations **myAlgEq, timeVariantAssignments **timeVarAssign, double sim_time, double dt, double *time, allocated_memory *mem, copied_AST *cp_AST){
  int i, j, k, l;
  int flag;
  //num of each objects
  int num_of_species = Model_getNumSpecies(m);
  int num_of_parameters = Model_getNumParameters(m);
  int num_of_compartments = Model_getNumCompartments(m);
  int num_of_reactions = Model_getNumReactions(m);
  int num_of_rules = Model_getNumRules(m);
  int num_of_events = Model_getNumEvents(m);
  int num_of_initialAssignments = Model_getNumInitialAssignments(m);
  ASTNode_t *node;

  //create mySpecies  
  Species_t *sp;
  for(i=0; i<num_of_species; i++){
    sp = (Species_t*)ListOf_get(Model_getListOfSpecies(m), i);
    mySp[i] = (mySpecies*)malloc(sizeof(mySpecies));
    mySp[i]->origin = sp;
    if(Species_isSetInitialAmount(sp)){
      mySp[i]->value = Species_getInitialAmount(sp);
      mySp[i]->is_amount = true;
      mySp[i]->is_concentration = false;
    }else if(Species_isSetInitialConcentration(sp)){
      mySp[i]->value = Species_getInitialConcentration(sp);
      mySp[i]->is_amount = false;
      mySp[i]->is_concentration = true;
    }else if(Species_getHasOnlySubstanceUnits(sp) || Compartment_getSpatialDimensions(Model_getCompartmentById(m, Species_getCompartment(sp))) == 0){
      mySp[i]->value = 0;
      mySp[i]->is_amount = true;
      mySp[i]->is_concentration = false;      
    }else{
      mySp[i]->value = 0;
      mySp[i]->is_amount = false;
      mySp[i]->is_concentration = true;
    }
    mySp[i]->has_only_substance_units = Species_getHasOnlySubstanceUnits(sp);
    mySp[i]->temp_value = mySp[i]->value;
    mySp[i]->locating_compartment = NULL;
    mySp[i]->delay_val = NULL;
    mySp[i]->depending_rule = NULL;
    mySp[i]->k[0] = 0;
    mySp[i]->k[1] = 0;
    mySp[i]->k[2] = 0;
    mySp[i]->k[3] = 0;
    mySp[i]->prev_val[0] = mySp[i]->value;
    mySp[i]->prev_val[1] = mySp[i]->value;
    mySp[i]->prev_val[2] = mySp[i]->value;
    mySp[i]->prev_k[0] = 0;
    mySp[i]->prev_k[1] = 0;
    mySp[i]->prev_k[2] = 0;
  }

  //create myParameters
  Parameter_t *param;
  for(i=0; i<num_of_parameters; i++){
    param = (Parameter_t*)ListOf_get(Model_getListOfParameters(m), i);
    myParam[i] = (myParameter*)malloc(sizeof(myParameter));
    myParam[i]->origin = param;
    if(Parameter_isSetValue(param)){
      myParam[i]->value = Parameter_getValue(param);
    }else{
      myParam[i]->value = 0;
    }
    myParam[i]->temp_value = myParam[i]->value;
    myParam[i]->delay_val = NULL;
    myParam[i]->depending_rule = NULL;
    myParam[i]->k[0] = 0;
    myParam[i]->k[1] = 0;
    myParam[i]->k[2] = 0;
    myParam[i]->k[3] = 0;
    myParam[i]->prev_val[0] = myParam[i]->value;
    myParam[i]->prev_val[1] = myParam[i]->value;
    myParam[i]->prev_val[2] = myParam[i]->value;
    myParam[i]->prev_k[0] = 0;
    myParam[i]->prev_k[1] = 0;
    myParam[i]->prev_k[2] = 0;
  }

  //create myCompartment
  Compartment_t *comp;
  for(i=0; i<num_of_compartments; i++){
    comp = (Compartment_t*)ListOf_get(Model_getListOfCompartments(m), i);
    myComp[i] = (myCompartment*)malloc(sizeof(myCompartment));
    myComp[i]->origin = comp;
    if(Compartment_isSetSize(comp)){
      myComp[i]->value = Compartment_getSize(comp);
    }else{
      myComp[i]->value = 1.0;
    }
    myComp[i]->temp_value = myComp[i]->value;
    myComp[i]->delay_val = NULL;
    myComp[i]->depending_rule = NULL;
    myComp[i]->k[0] = 0;
    myComp[i]->k[1] = 0;
    myComp[i]->k[2] = 0;
    myComp[i]->k[3] = 0;
    myComp[i]->prev_val[0] = myComp[i]->value;
    myComp[i]->prev_val[1] = myComp[i]->value;
    myComp[i]->prev_val[2] = myComp[i]->value;
    myComp[i]->prev_k[0] = 0;
    myComp[i]->prev_k[1] = 0;
    myComp[i]->prev_k[2] = 0;
    myComp[i]->num_of_including_species = 0;
  }

  //determine species locating compartment
  for(i=0; i<num_of_species; i++){
    for(j=0; j<num_of_compartments; j++){
      if(strcmp(Species_getCompartment(mySp[i]->origin), Compartment_getId(myComp[j]->origin)) == 0){
        mySp[i]->locating_compartment = myComp[j];
        myComp[j]->including_species[myComp[j]->num_of_including_species] = mySp[i];
        myComp[j]->num_of_including_species++;
      }
    }
  }

  //create myReaction & mySpeciseReference without equation
  Reaction_t *re;
  const char *product_id, *reactant_id;
  for(i=0; i<num_of_reactions; i++){
    re = (Reaction_t*)ListOf_get(Model_getListOfReactions(m), i);
    myRe[i] = (myReaction*)malloc(sizeof(myReaction));
    myRe[i]->origin = re;
    myRe[i]->is_fast = Reaction_getFast(re);
    myRe[i]->is_reversible = Reaction_getReversible(re);
    myRe[i]->num_of_products = 0;
    myRe[i]->num_of_reactants = 0;
    myRe[i]->products = (mySpeciesReference**)malloc(sizeof(mySpeciesReference*)*Reaction_getNumProducts(myRe[i]->origin));
    myRe[i]->reactants = (mySpeciesReference**)malloc(sizeof(mySpeciesReference*)*Reaction_getNumReactants(myRe[i]->origin));
    //products start
    for(j=0; j<Reaction_getNumProducts(myRe[i]->origin); j++){
      product_id = SpeciesReference_getSpecies((SpeciesReference_t*)ListOf_get(Reaction_getListOfProducts(myRe[i]->origin), j));
      for(k=0; k<num_of_species; k++){
        if(strcmp(product_id, Species_getId(mySp[k]->origin)) == 0){
          myRe[i]->products[myRe[i]->num_of_products] = (mySpeciesReference*)malloc(sizeof(mySpeciesReference));
          myRe[i]->products[myRe[i]->num_of_products]->origin = (SpeciesReference_t*)ListOf_get(Reaction_getListOfProducts(myRe[i]->origin), j);
          myRe[i]->products[myRe[i]->num_of_products]->mySp = mySp[k];
          myRe[i]->products[myRe[i]->num_of_products]->delay_val = NULL;
          myRe[i]->products[myRe[i]->num_of_products]->depending_rule = NULL;
          myRe[i]->num_of_products++;
        }
      }
    }//products fin
    //reactants start
    for(j=0; j<Reaction_getNumReactants(myRe[i]->origin); j++){
      reactant_id = SpeciesReference_getSpecies((SpeciesReference_t*)ListOf_get(Reaction_getListOfReactants(myRe[i]->origin), j));
      for(k=0; k<num_of_species; k++){
        if(strcmp(reactant_id, Species_getId(mySp[k]->origin)) == 0){
          myRe[i]->reactants[myRe[i]->num_of_reactants] = (mySpeciesReference*)malloc(sizeof(mySpeciesReference));
          myRe[i]->reactants[myRe[i]->num_of_reactants]->origin = (SpeciesReference_t*)ListOf_get(Reaction_getListOfReactants(myRe[i]->origin), j);
          myRe[i]->reactants[myRe[i]->num_of_reactants]->mySp = mySp[k];
          myRe[i]->reactants[myRe[i]->num_of_reactants]->delay_val = NULL;
          myRe[i]->reactants[myRe[i]->num_of_reactants]->depending_rule = NULL;
          myRe[i]->num_of_reactants++;
        }
      }
    }//reactants fin
    myRe[i]->products_equili_numerator = NULL;
    myRe[i]->reactants_equili_numerator = NULL;
  }

  //create myInitialAssignments
  InitialAssignment_t *initAssign;
  char *time_variant_target_id[MAX_DELAY_REACTION_NUM];
  int num_of_time_variant_targets = 0;
  for(i=0; i<num_of_initialAssignments; i++){
    initAssign = (InitialAssignment_t*)ListOf_get(Model_getListOfInitialAssignments(m), i);
    myInitAssign[i] = (myInitialAssignment*)malloc(sizeof(myInitialAssignment));
    myInitAssign[i]->origin = initAssign;
    myInitAssign[i]->target_species = NULL;
    myInitAssign[i]->target_parameter = NULL;
    myInitAssign[i]->target_compartment = NULL;
    myInitAssign[i]->target_species_reference = NULL;
    flag = 1;
    for(j=0; j<num_of_species; j++){
      if(strcmp(InitialAssignment_getSymbol(initAssign), Species_getId(mySp[j]->origin)) == 0){
        myInitAssign[i]->target_species = mySp[j];
        flag = 0;
        break;
      }
    }
    if(flag){
      for(j=0; j<num_of_parameters; j++){
        if(strcmp(InitialAssignment_getSymbol(initAssign), Parameter_getId(myParam[j]->origin)) == 0){
          myInitAssign[i]->target_parameter = myParam[j];
          flag = 0;
          break;
        }
      }
    }
    if(flag){
      for(j=0; j<num_of_compartments; j++){
        if(strcmp(InitialAssignment_getSymbol(initAssign), Compartment_getId(myComp[j]->origin)) == 0){
          myInitAssign[i]->target_compartment = myComp[j];
          flag = 0;
          break;
        }
      }
    }
    node = (ASTNode_t*)InitialAssignment_getMath(initAssign);
    node = ASTNode_deepCopy(node);
    dbg_printf("original math : ");
    check_AST(node, NULL);
    //alter_tree_structure(m, &node, cp_AST);
    alter_tree_structure(m, &node, NULL, 0, cp_AST);
    //unit
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
    //unit
    dbg_printf("altered math : ");
    check_AST(node, NULL);
    if(include_time(node, 0)){
      time_variant_target_id[num_of_time_variant_targets++] = (char*)InitialAssignment_getSymbol(initAssign);
    }
    myInitAssign[i]->eq = (equation*)malloc(sizeof(equation));
    myInitAssign[i]->eq->math_length = get_equation(m, myInitAssign[i]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, NULL, time_variant_target_id, 0, NULL, mem);
    //ASTNode_free(node);
    dbg_printf("math\n");
    check_math(myInitAssign[i]->eq);
  }

  //find time variant target of assignment rule
  Rule_t *rule;
  *timeVarAssign = (timeVariantAssignments*)malloc(sizeof(timeVariantAssignments));
  (*timeVarAssign)->num_of_time_variant_assignments = 0;
  for(i=0; i<num_of_rules; i++){
    rule = Model_getRule(m, i);
    node = (ASTNode_t*)Rule_getMath(rule);
    if(Rule_isAssignment(rule) && include_time(node, 0)){
      node = ASTNode_deepCopy(node);
      (*timeVarAssign)->target_id[(*timeVarAssign)->num_of_time_variant_assignments] = (char*)Rule_getVariable(rule);
      dbg_printf("original math : ");
      check_AST(node, NULL);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      //unit
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
      //unit
      dbg_printf("altered math : ");
      check_AST(node, NULL);
      (*timeVarAssign)->eq[(*timeVarAssign)->num_of_time_variant_assignments] = (equation*)malloc(sizeof(equation));
      (*timeVarAssign)->eq[(*timeVarAssign)->num_of_time_variant_assignments]->math_length = get_equation(m, (*timeVarAssign)->eq[(*timeVarAssign)->num_of_time_variant_assignments], mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, NULL, time_variant_target_id, 0, NULL, mem);
      ((*timeVarAssign)->num_of_time_variant_assignments)++;
    }
  }

  //create myReaction & mySpeciseReference equation
  ASTNode_t *times_node, *conv_factor_node;
  for(i=0; i<num_of_reactions; i++){
    re = (Reaction_t*)ListOf_get(Model_getListOfReactions(m), i);
    node = (ASTNode_t*)KineticLaw_getMath(Reaction_getKineticLaw(re));
    node = ASTNode_deepCopy(node);
    dbg_printf("original math of %s: ", Reaction_getId(myRe[i]->origin));
    check_AST(node, NULL);
    alter_tree_structure(m, &node, NULL, 0, cp_AST);
    set_local_para_as_value(node, Reaction_getKineticLaw(re));
    dbg_printf("alterated math of %s : ", Reaction_getId(myRe[i]->origin));
    check_AST(node, NULL);
    myRe[i]->eq = (equation*)malloc(sizeof(equation));
    myRe[i]->eq->math_length = get_equation(m, myRe[i]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
    //ASTNode_free(node);
    dbg_printf("math of %s\n", Reaction_getId(myRe[i]->origin));
    check_math(myRe[i]->eq);
    //products start
    for(j=0; j<myRe[i]->num_of_products; j++){
      myRe[i]->products[j]->eq = (equation*)malloc(sizeof(equation));
      myRe[i]->products[j]->eq->math_length = 0;
      node = (ASTNode_t*)StoichiometryMath_getMath(SpeciesReference_getStoichiometryMath(myRe[i]->products[j]->origin));
      if(node != NULL){//l2v4
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
        if(isnan(SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin)) && SpeciesReference_isSetId(myRe[i]->products[j]->origin)){//l3v1
          node = ASTNode_createWithType(AST_NAME);
          ASTNode_setName(node, SpeciesReference_getId(myRe[i]->products[j]->origin));
          myRe[i]->products[j]->value = 1.0;
          myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
        }else if(isnan(SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin))){//l3v1
          node = ASTNode_createWithType(AST_REAL);
          ASTNode_setReal(node, 1.0);
          myRe[i]->products[j]->value = 1.0;
          myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
        }else{//l2v4
          node = ASTNode_createWithType(AST_REAL);
          ASTNode_setReal(node, SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin));
          myRe[i]->products[j]->value = SpeciesReference_getStoichiometry(myRe[i]->products[j]->origin);
          myRe[i]->products[j]->temp_value = myRe[i]->products[j]->value;
        }
      }
      dbg_printf("original math of %s : ", SpeciesReference_getSpecies(myRe[i]->products[j]->origin));
      check_AST(node, NULL);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      //unit
      if(myRe[i]->products[j]->mySp->is_concentration){
        assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myRe[i]->products[j]->mySp->origin))), 1);
      }
      //unit
      //conversion factor
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
      //conversion factor
      myRe[i]->products[j]->eq->math_length = get_equation(m, myRe[i]->products[j]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
      dbg_printf("alterated math of %s : ", SpeciesReference_getSpecies(myRe[i]->products[j]->origin));
      check_AST(node, NULL);
      //ASTNode_free(node);
      dbg_printf("math of %s\n", SpeciesReference_getSpecies(myRe[i]->products[j]->origin));
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
    //products fin
    //reactants start
    for(j=0; j<myRe[i]->num_of_reactants; j++){
      myRe[i]->reactants[j]->depending_rule = NULL;
      myRe[i]->reactants[j]->eq = (equation*)malloc(sizeof(equation));
      myRe[i]->reactants[j]->eq->math_length = 0;
      node = (ASTNode_t*)StoichiometryMath_getMath(SpeciesReference_getStoichiometryMath(myRe[i]->reactants[j]->origin));
      if(node != NULL){//l2v4
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
        if(isnan(SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin)) && SpeciesReference_isSetId(myRe[i]->reactants[j]->origin)){//l3v1
          node = ASTNode_createWithType(AST_NAME);
          ASTNode_setName(node, SpeciesReference_getId(myRe[i]->reactants[j]->origin));
          myRe[i]->reactants[j]->value = 1.0;
          myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
        }else if(isnan(SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin))){//l3v1
          node = ASTNode_createWithType(AST_REAL);
          ASTNode_setReal(node, 1.0);
          myRe[i]->reactants[j]->value = 1.0;
          myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
        }else{//l2v4
          node = ASTNode_createWithType(AST_REAL);
          ASTNode_setReal(node, SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin));
          myRe[i]->reactants[j]->value = SpeciesReference_getStoichiometry(myRe[i]->reactants[j]->origin);
          myRe[i]->reactants[j]->temp_value = myRe[i]->reactants[j]->value;
        }
      }
      dbg_printf("original math of %s : ", SpeciesReference_getSpecies(myRe[i]->reactants[j]->origin));
      check_AST(node, NULL);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      //unit
      if(myRe[i]->reactants[j]->mySp->is_concentration){
        assignment_alter_tree_structure(&node, (char*)Compartment_getId(Model_getCompartmentById(m, Species_getCompartment(myRe[i]->reactants[j]->mySp->origin))), 1);
      }
      //unit
      //conversion factor
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
      //conversion factor
      myRe[i]->reactants[j]->eq->math_length = get_equation(m, myRe[i]->reactants[j]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
      dbg_printf("altered math of %s : ", SpeciesReference_getSpecies(myRe[i]->reactants[j]->origin));
      check_AST(node, NULL);
      //ASTNode_free(node);
      dbg_printf("math of %s\n", SpeciesReference_getSpecies(myRe[i]->reactants[j]->origin));
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
      //reactants fin
    }
  }


  //Initial Assignment target_species_reference
  for(i=0; i<num_of_initialAssignments; i++){
    flag = 1;
    if(flag){
      for(j=0; j<num_of_reactions; j++){
        for(k=0; k<myRe[j]->num_of_products; k++){
          if(SpeciesReference_isSetId(myRe[j]->products[k]->origin) 
              && strcmp(InitialAssignment_getSymbol(myInitAssign[i]->origin), SpeciesReference_getId(myRe[j]->products[k]->origin)) == 0){
            myInitAssign[i]->target_species_reference = myRe[j]->products[k];
            flag = 0;
            break;
          }
        }
        if(!flag){
          break;
        }
        for(k=0; k<myRe[j]->num_of_reactants; k++){
          if(SpeciesReference_isSetId(myRe[j]->reactants[k]->origin) 
              && strcmp(InitialAssignment_getSymbol(myInitAssign[i]->origin), SpeciesReference_getId(myRe[j]->reactants[k]->origin)) == 0){
            myInitAssign[i]->target_species_reference = myRe[j]->reactants[k];
            flag = 0;
            break;
          }
        }
        if(!flag){
          break;
        }
      }
    }
  }

  //prepare reversible fast reaction
  prepare_reversible_fast_reaction(m, myRe, mySp, myParam, myComp, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem, cp_AST);

  //create myRules
  for(i=0; i<num_of_rules; i++){
    rule = (Rule_t*)ListOf_get(Model_getListOfRules(m), i);
    myRu[i] = (myRule*)malloc(sizeof(myRule));
    myRu[i]->origin = rule;
    myRu[i]->target_species = NULL;
    myRu[i]->target_parameter = NULL;
    myRu[i]->target_compartment = NULL;
    myRu[i]->target_species_reference = NULL;
    myRu[i]->is_rate = Rule_isRate(rule);
    myRu[i]->is_assignment = Rule_isAssignment(rule);
    myRu[i]->is_algebraic = Rule_isAlgebraic(rule);
    if(Rule_isRate(rule) || Rule_isAssignment(rule)){
      flag = 1;
      for(j=0; j<num_of_species; j++){
        if(strcmp(Rule_getVariable(rule), Species_getId(mySp[j]->origin)) == 0){
          myRu[i]->target_species = mySp[j];
          mySp[j]->depending_rule = myRu[i];
          flag = 0;
          break;
        }
      }
      if(flag){
        for(j=0; j<num_of_parameters; j++){
          if(strcmp(Rule_getVariable(rule), Parameter_getId(myParam[j]->origin)) == 0){
            myRu[i]->target_parameter = myParam[j];
            myParam[j]->depending_rule = myRu[i];
            flag = 0;
            break;
          }
        }
      }
      if(flag){
        for(j=0; j<num_of_compartments; j++){
          if(strcmp(Rule_getVariable(rule), Compartment_getId(myComp[j]->origin)) == 0){
            myRu[i]->target_compartment = myComp[j];
            myComp[j]->depending_rule = myRu[i];
            flag = 0;
            break;
          }
        }
      }
      if(flag){
        for(j=0; j<num_of_reactions; j++){
          for(k=0; k<myRe[j]->num_of_products; k++){
            if(SpeciesReference_isSetId(myRe[j]->products[k]->origin) 
                && strcmp(Rule_getVariable(rule), SpeciesReference_getId(myRe[j]->products[k]->origin)) == 0){
              myRu[i]->target_species_reference = myRe[j]->products[k];
              myRe[j]->products[k]->depending_rule = myRu[i];
              flag = 0;
              break;
            }
          }
          if(!flag){
            break;
          }
          for(k=0; k<myRe[j]->num_of_reactants; k++){
            if(SpeciesReference_isSetId(myRe[j]->reactants[k]->origin) 
                && strcmp(Rule_getVariable(rule), SpeciesReference_getId(myRe[j]->reactants[k]->origin)) == 0){
              myRu[i]->target_species_reference = myRe[j]->reactants[k];
              myRe[j]->reactants[k]->depending_rule = myRu[i];
              flag = 0;
              break;
            }
          }
          if(!flag){
            break;
          }
        }
      }
      node = (ASTNode_t*)Rule_getMath(rule);
      node = ASTNode_deepCopy(node);
      dbg_printf("original math : ");
      check_AST(node, NULL);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      //unit
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
      //unit
      dbg_printf("altered math : ");
      check_AST(node, NULL);
      myRu[i]->eq = (equation*)malloc(sizeof(equation));
      myRu[i]->eq->math_length = get_equation(m, myRu[i]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
      //ASTNode_free(node);
      dbg_printf("math\n");
      check_math(myRu[i]->eq);
    }
  }

  //create myEvents & myEventAssignments
  Event_t *event;
  for(i=0; i<num_of_events; i++){
    event = (Event_t*)ListOf_get(Model_getListOfEvents(m), i);
    myEv[i] = (myEvent*)malloc(sizeof(myEvent));
    myEv[i]->origin = event;
    if(Trigger_isSetPersistent(Event_getTrigger(event))){
      myEv[i]->is_persistent = Trigger_getPersistent(Event_getTrigger(event));
    }else{
      myEv[i]->is_persistent = true;
    }
    node = (ASTNode_t*)Trigger_getMath(Event_getTrigger(event));
    node = ASTNode_deepCopy(node);
    dbg_printf("original math of %s : ", Event_getId(myEv[i]->origin));
    check_AST(node, NULL);
    //alter_tree_structure(m, &node, cp_AST);
    alter_tree_structure(m, &node, NULL, 0, cp_AST);
    myEv[i]->eq = (equation*)malloc(sizeof(equation));
    myEv[i]->eq->math_length = get_equation(m, myEv[i]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
    dbg_printf("altered math of %s : ", Event_getId(myEv[i]->origin));
    check_AST(node, NULL);
    //ASTNode_free(node);
    dbg_printf("math of %s\n", Event_getId(myEv[i]->origin));
    check_math(myEv[i]->eq);
    if(Trigger_isSetInitialValue(Event_getTrigger(event))){
      myEv[i]->is_able_to_fire = !Trigger_getInitialValue(Event_getTrigger(event));
    }else{
      myEv[i]->is_able_to_fire = true;
    }
    myEv[i]->assignments = (myEventAssignment**)malloc(sizeof(myEventAssignment*)*Event_getNumEventAssignments(event));
    for(j=0; j<Event_getNumEventAssignments(event); j++){
      myEv[i]->assignments[j] = (myEventAssignment*)malloc(sizeof(myEventAssignment));
      myEv[i]->assignments[j]->origin = (EventAssignment_t*)ListOf_get(Event_getListOfEventAssignments(event), j);
      myEv[i]->assignments[j]->target_species = NULL;
      myEv[i]->assignments[j]->target_parameter = NULL;
      myEv[i]->assignments[j]->target_compartment = NULL;
      myEv[i]->assignments[j]->target_species_reference = NULL;
      flag = 1;
      for(k=0; k<num_of_species; k++){
        if(strcmp(EventAssignment_getVariable(myEv[i]->assignments[j]->origin), Species_getId(mySp[k]->origin)) == 0){
          myEv[i]->assignments[j]->target_species = mySp[k];
          flag = 0;
          break;
        }
      }
      if(flag){
        for(k=0; k<num_of_parameters; k++){
          if(strcmp(EventAssignment_getVariable(myEv[i]->assignments[j]->origin), Parameter_getId(myParam[k]->origin)) == 0){
            myEv[i]->assignments[j]->target_parameter = myParam[k];
            flag = 0;
            break;
          }
        }
      }
      if(flag){
        for(k=0; k<num_of_compartments; k++){
          if(strcmp(EventAssignment_getVariable(myEv[i]->assignments[j]->origin), Compartment_getId(myComp[k]->origin)) == 0){
            myEv[i]->assignments[j]->target_compartment = myComp[k];
            flag = 0;
            break;
          }
        }
      }
      if(flag){
        for(k=0; k<num_of_reactions; k++){
          for(l=0; l<myRe[k]->num_of_products; l++){
            if(SpeciesReference_isSetId(myRe[k]->products[l]->origin) 
                && strcmp(EventAssignment_getVariable(myEv[i]->assignments[j]->origin), SpeciesReference_getId(myRe[k]->products[l]->origin)) == 0){
              myEv[i]->assignments[j]->target_species_reference = myRe[k]->products[l];
              flag = 0;
              break;
            }
          }
          if(!flag){
            break;
          }
          for(l=0; l<myRe[k]->num_of_reactants; l++){
            if(SpeciesReference_isSetId(myRe[k]->reactants[l]->origin) 
                && strcmp(EventAssignment_getVariable(myEv[i]->assignments[j]->origin), SpeciesReference_getId(myRe[k]->reactants[l]->origin)) == 0){
              myEv[i]->assignments[j]->target_species_reference = myRe[k]->reactants[l];
              flag = 0;
              break;
            }
          }
          if(!flag){
            break;
          }
        }
      }
      node = (ASTNode_t*)EventAssignment_getMath(myEv[i]->assignments[j]->origin);
      node = ASTNode_deepCopy(node);
      dbg_printf("original math : ");
      check_AST(node, NULL);
      if(Event_getDelay(myEv[i]->origin) != NULL && Event_getUseValuesFromTriggerTime(myEv[i]->origin)){
        pre_ev_alter_tree_structure(&node, NULL, 0, (ASTNode_t*)Delay_getMath(Event_getDelay(myEv[i]->origin)));
        dbg_printf("after pre ev alter math : ");
        check_AST(node, NULL);
        ev_alter_tree_structure(m, &node, NULL, 0, cp_AST);
        dbg_printf("after ev alter math : ");
        check_AST(node, NULL);
        post_ev_alter_tree_structure(m, &node, NULL, 0);
        dbg_printf("after post ev alter math : ");
        check_AST(node, NULL);
      }else{
        alter_tree_structure(m, &node, NULL, 0, cp_AST);
      }
      //unit
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
      //unit
      dbg_printf("altered math : ");
      check_AST(node, NULL);
      myEv[i]->assignments[j]->eq = (equation*)malloc(sizeof(equation));
      myEv[i]->assignments[j]->eq->math_length = get_equation(m, myEv[i]->assignments[j]->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
      //ASTNode_free(node);
      dbg_printf("math\n");
      check_math(myEv[i]->assignments[j]->eq);
    }
    myEv[i]->event_delay = NULL;
    if(Event_getDelay(event) != NULL){
      myEv[i]->event_delay = (myDelay*)malloc(sizeof(myDelay));
      myEv[i]->event_delay->origin = Event_getDelay(event);
      node = (ASTNode_t*)Delay_getMath(myEv[i]->event_delay->origin);
      node = ASTNode_deepCopy(node);
      dbg_printf("original math : ");
      check_AST(node, NULL);
      //alter_tree_structure(m, &node, cp_AST);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      myEv[i]->event_delay->eq = (equation*)malloc(sizeof(equation));
      myEv[i]->event_delay->eq->math_length = get_equation(m, myEv[i]->event_delay->eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
      dbg_printf("altered math : ");
      check_AST(node, NULL);
      //ASTNode_free(node);
      dbg_printf("math\n");
      check_math(myEv[i]->event_delay->eq);
      myEv[i]->firing_times = (double*)malloc(sizeof(double)*(int)(sim_time/dt));
      for(j=0; j<(int)(sim_time/dt); j++){
        myEv[i]->firing_times[j] = sim_time+1;
      }
      myEv[i]->num_of_delayed_events_que = 0;
      myEv[i]->next_firing_index = 0;
    }
    myEv[i]->priority_eq = NULL;
    if(Event_isSetPriority(myEv[i]->origin)){
      myEv[i]->priority_eq = (equation*)malloc(sizeof(equation));
      node = (ASTNode_t*)Priority_getMath(Event_getPriority(myEv[i]->origin));
      node = ASTNode_deepCopy(node);
      alter_tree_structure(m, &node, NULL, 0, cp_AST);
      myEv[i]->priority_eq->math_length = get_equation(m, myEv[i]->priority_eq, mySp, myParam, myComp, myRe, node, 0, sim_time, dt, time, myInitAssign, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem);
    }
  }

  //prepare Algebraic
  int num_of_algebraic_rules = 0;
  myAlgebraicEquations *algEq;
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
          algEq->coefficient_matrix[i][j] = (equation*)malloc(sizeof(equation));
          algEq->coefficient_matrix[i][j]->math_length = 0;
        }
      }
      algEq->constant_vector = (equation**)malloc(sizeof(equation*)*algEq->num_of_algebraic_rules);
      for(i=0; i<algEq->num_of_algebraic_rules; i++){
        algEq->constant_vector[i] = (equation*)malloc(sizeof(equation));
      }
    }else{
      algEq->coefficient = (equation*)malloc(sizeof(equation));
      algEq->constant = (equation*)malloc(sizeof(equation));
    }
    dbg_printf("prepare algebraic start\n");
    prepare_algebraic(m, mySp, myParam, myComp, myRe, myRu, myEv, myInitAssign, algEq, sim_time, dt, time, time_variant_target_id, num_of_time_variant_targets, *timeVarAssign, mem, cp_AST);
    dbg_printf("prepare algebraic finish\n");
    if(algEq->num_of_algebraic_rules > 1){
      for(i=0; i<algEq->num_of_algebraic_rules; i++){
        for(j=0; j<algEq->num_of_algebraic_rules; j++){
          dbg_printf("math of coefficient matrix[%d][%d] is\n", i, j);
          check_math(algEq->coefficient_matrix[i][j]);
        }
      }
      for(i=0; i<algEq->num_of_algebraic_rules; i++){
        dbg_printf("constant vector[%d] is\n", i);
        check_math(algEq->constant_vector[i]);
      }
    }else{
      dbg_printf("math of coefficient is\n");
      check_math(algEq->coefficient);
      dbg_printf("math of constant is\n");
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
