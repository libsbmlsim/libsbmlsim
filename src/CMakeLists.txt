set (SOURCES_COMMON)
set(SOURCES_COMMON ${SOURCES_COMMON}
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/common.h
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/libsbmlsim.h
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/methods.h
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/myResult.h
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/my_getopt.h
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/osarch.h
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/version.h
  ${PROJECT_SOURCE_DIR}/src/alg_alter_tree_structure.c
  ${PROJECT_SOURCE_DIR}/src/alter_tree_structure.c
  ${PROJECT_SOURCE_DIR}/src/assignment_alter_tree_structure.c
  ${PROJECT_SOURCE_DIR}/src/calc.c
  ${PROJECT_SOURCE_DIR}/src/calc_event.c
  ${PROJECT_SOURCE_DIR}/src/calc_initial_assignment.c
  ${PROJECT_SOURCE_DIR}/src/calc_k.c
  ${PROJECT_SOURCE_DIR}/src/calc_temp_value.c
  ${PROJECT_SOURCE_DIR}/src/check_AST.c
  ${PROJECT_SOURCE_DIR}/src/check_math.c
  ${PROJECT_SOURCE_DIR}/src/check_num.c
  ${PROJECT_SOURCE_DIR}/src/chomp.c
  ${PROJECT_SOURCE_DIR}/src/create_calc_object_list.c
  ${PROJECT_SOURCE_DIR}/src/create_myResult_content.c
  ${PROJECT_SOURCE_DIR}/src/create_mySBML_objects.c
  ${PROJECT_SOURCE_DIR}/src/dbg_printf.c
  ${PROJECT_SOURCE_DIR}/src/e_acosh.c
  ${PROJECT_SOURCE_DIR}/src/e_atanh.c
  ${PROJECT_SOURCE_DIR}/src/ev_alter_tree_structure.c
  ${PROJECT_SOURCE_DIR}/src/forwarding_value.c
  ${PROJECT_SOURCE_DIR}/src/free_mySBML_objects.c
  ${PROJECT_SOURCE_DIR}/src/get_equation.c
  ${PROJECT_SOURCE_DIR}/src/initialize_delay_val.c
  ${PROJECT_SOURCE_DIR}/src/lu_decomposition.c
  ${PROJECT_SOURCE_DIR}/src/lu_solve.c
  ${PROJECT_SOURCE_DIR}/src/math_functions.c
  ${PROJECT_SOURCE_DIR}/src/myASTNode_func.c
  ${PROJECT_SOURCE_DIR}/src/output_result.c
  ${PROJECT_SOURCE_DIR}/src/prepare_algebraic.c
  ${PROJECT_SOURCE_DIR}/src/prepare_reversible_fast_reaction.c
  ${PROJECT_SOURCE_DIR}/src/prg_printf.c
  ${PROJECT_SOURCE_DIR}/src/print_node_type.c
  ${PROJECT_SOURCE_DIR}/src/print_result_list.c
  ${PROJECT_SOURCE_DIR}/src/s_log1p.c
  ${PROJECT_SOURCE_DIR}/src/set_local_para_as_value.c
  ${PROJECT_SOURCE_DIR}/src/simulate_explicit.c
  ${PROJECT_SOURCE_DIR}/src/simulate_implicit.c
  ${PROJECT_SOURCE_DIR}/src/str_is_number.c
  ${PROJECT_SOURCE_DIR}/src/substitute_delay_val.c
  ${PROJECT_SOURCE_DIR}/src/util.c
  )

set (SOURCES_LIB ${SOURCES_COMMON} ${PROJECT_SOURCE_DIR}/src/lib_main.c)
set (SOURCES_MAIN ${SOURCES_COMMON} ${PROJECT_SOURCE_DIR}/src/main.c)
set (SOURCES_TEST test.c)

include_directories ("${PROJECT_SOURCE_DIR}/src/libsbmlsim" ${LIBSBML_INCLUDE_DIR})

# Static library (only for UNIX)
#if(NOT MSVC)
  add_library(sbmlsim-static STATIC ${SOURCES_LIB})
  target_link_libraries(sbmlsim-static ${LIBSBML_LIBRARIES})
#endif()

# Shared library
add_library(sbmlsim SHARED ${SOURCES_LIB})
target_link_libraries(sbmlsim ${LIBSBML_LIBRARIES})

# Test program
add_executable(test ${SOURCES_TEST})
target_link_libraries(test sbmlsim ${LIBSBML_LIBRARIES})

# Simulator for test suiteprogram
add_executable(simulateSBML ${SOURCES_MAIN})
target_link_libraries(simulateSBML sbmlsim-static ${LIBSBML_LIBRARIES})
if(NOT MSVC)
ADD_CUSTOM_COMMAND(TARGET simulateSBML
  POST_BUILD
  COMMAND cmake -E echo "Copying simulator to testcase directory..."
  COMMAND cmake -E copy simulateSBML ${PROJECT_SOURCE_DIR}/testcases
  )
endif()

# Install targets
INSTALL(TARGETS simulateSBML sbmlsim-static sbmlsim
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )

# Install include files
FILE(GLOB include_files "${PROJECT_SOURCE_DIR}/src/libsbmlsim/*.h")
INSTALL(FILES ${include_files} DESTINATION include/libsbmlsim)
# Install documents
FILE(GLOB text_files "${PROJECT_SOURCE_DIR}/*.txt")
INSTALL(FILES ${text_files} DESTINATION share/libsbmlsim)
# Install sample SBML, result files
FILE(GLOB sample_files "${PROJECT_SOURCE_DIR}/examples/sample.*")
INSTALL(FILES ${sample_files} DESTINATION share/libsbmlsim)
# Install C and C++ code
INSTALL(FILES ${PROJECT_SOURCE_DIR}/examples/c/test.c ${PROJECT_SOURCE_DIR}/examples/c/Makefile DESTINATION share/libsbmlsim/c)
INSTALL(FILES ${PROJECT_SOURCE_DIR}/examples/cpp/test.cpp ${PROJECT_SOURCE_DIR}/examples/cpp/Makefile DESTINATION share/libsbmlsim/cpp)

# Bindings
if(WITH_JAVA)
add_subdirectory(bindings/java)
endif()

if(WITH_PYTHON)
add_subdirectory(bindings/python)
endif()

if(WITH_PERL)
add_subdirectory(bindings/perl)
endif()

if(WITH_RUBY)
add_subdirectory(bindings/ruby)
endif()


