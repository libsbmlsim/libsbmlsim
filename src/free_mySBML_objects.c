#include "sbml/common/common.h"
#include "sbml/SBMLReader.h"
#include "sbml/SBMLTypes.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "header.h"

void free_mySBML_objects(Model_t *m, mySpecies *mySp[], myParameter *myParam[], myCompartment *myComp[], myReaction *myRe[], myRule *myRu[], myEvent *myEv[], myInitialAssignment *myInitAssign[], myAlgebraicEquations *myAlgEq, timeVariantAssignments *timeVarAssign, double sim_time, double dt, allocated_memory *mem, copied_AST *cp_AST){
  int i, j;

  //free
  for(i=0; i<Model_getNumSpecies(m); i++){
    if(mySp[i]->delay_val != NULL){
      for(j=0; j<(int)(sim_time/dt+1); j++){
	free(mySp[i]->delay_val[j]);
      }
      free(mySp[i]->delay_val);
    }
    free(mySp[i]);
  }
  for(i=0; i<Model_getNumParameters(m); i++){
    if(myParam[i]->delay_val != NULL){
      for(j=0; j<(int)(sim_time/dt+1); j++){
	free(myParam[i]->delay_val[j]);
      }
      free(myParam[i]->delay_val);
    }
    free(myParam[i]);
  }
  for(i=0; i<Model_getNumCompartments(m); i++){
    if(myComp[i]->delay_val != NULL){
      for(j=0; j<(int)(sim_time/dt+1); j++){
	free(myComp[i]->delay_val[j]);
      }
      free(myComp[i]->delay_val);
    }
    free(myComp[i]);
  }
  for(i=0; i<Model_getNumReactions(m); i++){
    for(j=0; j<myRe[i]->num_of_products; j++){
      free(myRe[i]->products[j]);
    }
    free(myRe[i]->products);
    for(j=0; j<myRe[i]->num_of_reactants; j++){
      free(myRe[i]->reactants[j]);
    }
    free(myRe[i]->reactants);
    free(myRe[i]->eq);
    if(myRe[i]->products_equili_numerator != NULL){
      free(myRe[i]->products_equili_numerator);
    }
    if(myRe[i]->reactants_equili_numerator != NULL){
      free(myRe[i]->reactants_equili_numerator);
    }
    free(myRe[i]);
  }
  for(i=0; i<Model_getNumRules(m); i++){
    if(!myRu[i]->is_algebraic){
      free(myRu[i]->eq);
    }
    free(myRu[i]);
  }
  for(i=0; i<Model_getNumEvents(m); i++){
    for(j=0; j<Event_getNumEventAssignments(myEv[i]->origin); j++){
      free(myEv[i]->assignments[j]->eq);
      free(myEv[i]->assignments[j]);
    }
    free(myEv[i]->assignments);
    if(myEv[i]->event_delay != NULL){
      free(myEv[i]->event_delay->eq);
      free(myEv[i]->event_delay);
      free(myEv[i]->firing_times);
    }
    if(myEv[i]->priority_eq != NULL){
      free(myEv[i]->priority_eq);
    }
    free(myEv[i]);
  }
  for(i=0; i<Model_getNumInitialAssignments(m); i++){
    free(myInitAssign[i]->eq);
    free(myInitAssign[i]);
  }
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
  for(i=0; i<mem->num_of_allocated_memory; i++){
    free(mem->memory[i]);
  }
  free(mem);
  //piecewiseを含むテストケースでsegmentation fault
/*   for(i=0; i<cp_AST->num_of_copied_AST; i++){ */
/*     ASTNode_free(cp_AST->ast[i]); */
/*   } */
/*   free(cp_AST); */
  printf("all allocated memory is freeed\n");
}
