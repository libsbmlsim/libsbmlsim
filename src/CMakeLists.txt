set (SOURCES_COMMON)
set(SOURCES_COMMON ${SOURCES_COMMON}
  libsbmlsim/common.h
  libsbmlsim/libsbmlsim.h
  libsbmlsim/methods.h
  libsbmlsim/version.h
  alg_alter_tree_structure.c
  alter_tree_structure.c
  assignment_alter_tree_structure.c
  calc.c
  calc_event.c
  calc_initial_assignment.c
  calc_k.c
  calc_temp_value.c
  check_AST.c
  check_math.c
  check_num.c
  chomp.c
  create_calc_object_list.c
  create_myResult_content.c
  create_mySBML_objects.c
  ev_alter_tree_structure.c
  forwarding_value.c
  free_mySBML_objects.c
  get_equation.c
  initialize_delay_val.c
  lu_decomposition.c
  lu_solve.c
  math_functions.c
  myASTNode_func.c
  output_result.c
  prepare_algebraic.c
  prepare_reversible_fast_reaction.c
  print_node_type.c
  print_result_list.c
  set_local_para_as_value.c
  simulate_explicit.c
  simulate_implicit.c
  str_is_number.c
  substitute_delay_val.c
  util.c
  )

set (SOURCES_LIB ${SOURCES_COMMON} lib_main.c)
set (SOURCES_MAIN ${SOURCES_COMMON} main.c)
set (SOURCES_TEST test.c)

include_directories("./libsbmlsim")

# Static library
add_library(sbmlsim-static STATIC ${SOURCES_LIB})
target_link_libraries(sbmlsim-static sbml)

# Shared library
add_library(sbmlsim SHARED ${SOURCES_LIB})
target_link_libraries(sbmlsim sbml)

# Test program
add_executable(test ${SOURCES_TEST})
target_link_libraries(test sbmlsim sbml)

