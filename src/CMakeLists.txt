set (SOURCES_COMMON)
set(SOURCES_COMMON ${SOURCES_COMMON}
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/common.h
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/libsbmlsim.h
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/dSFMT-common.h
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/dSFMT-params.h
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/dSFMT-params19937.h
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/dSFMT.h
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/methods.h
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/myResult.h
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/my_getopt.h
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/osarch.h
  ${PROJECT_SOURCE_DIR}/src/libsbmlsim/version.h
  ${PROJECT_SOURCE_DIR}/src/alg_alter_tree_structure.c
  ${PROJECT_SOURCE_DIR}/src/allocated_memory.c
  ${PROJECT_SOURCE_DIR}/src/alter_tree_structure.c
  ${PROJECT_SOURCE_DIR}/src/assignment_alter_tree_structure.c
  ${PROJECT_SOURCE_DIR}/src/bifurcation_analysis.c
  ${PROJECT_SOURCE_DIR}/src/check_AST.c
  ${PROJECT_SOURCE_DIR}/src/check_math.c
  ${PROJECT_SOURCE_DIR}/src/check_num.c
  ${PROJECT_SOURCE_DIR}/src/copied_AST.c
  ${PROJECT_SOURCE_DIR}/src/count_ode.c
  ${PROJECT_SOURCE_DIR}/src/dSFMT.c
  ${PROJECT_SOURCE_DIR}/src/ev_alter_tree_structure.c
  ${PROJECT_SOURCE_DIR}/src/equation.c
  ${PROJECT_SOURCE_DIR}/src/get_equation.c
  ${PROJECT_SOURCE_DIR}/src/myASTNode.c
  ${PROJECT_SOURCE_DIR}/src/myCompartment.c
  ${PROJECT_SOURCE_DIR}/src/myDelay.c
  ${PROJECT_SOURCE_DIR}/src/myEvent.c
  ${PROJECT_SOURCE_DIR}/src/myEventAssignment.c
  ${PROJECT_SOURCE_DIR}/src/myInitialAssignment.c
  ${PROJECT_SOURCE_DIR}/src/myParameter.c
  ${PROJECT_SOURCE_DIR}/src/myReaction.c
  ${PROJECT_SOURCE_DIR}/src/myResult.c
  ${PROJECT_SOURCE_DIR}/src/myRule.c
  ${PROJECT_SOURCE_DIR}/src/mySpecies.c
  ${PROJECT_SOURCE_DIR}/src/mySpeciesReference.c
  ${PROJECT_SOURCE_DIR}/src/mySBML_objects.c
  ${PROJECT_SOURCE_DIR}/src/output_result.c
  ${PROJECT_SOURCE_DIR}/src/prepare_algebraic.c
  ${PROJECT_SOURCE_DIR}/src/prepare_reversible_fast_reaction.c
  ${PROJECT_SOURCE_DIR}/src/print_node_type.c
  ${PROJECT_SOURCE_DIR}/src/print_result_list.c
  ${PROJECT_SOURCE_DIR}/src/search_max.c
  ${PROJECT_SOURCE_DIR}/src/set_local_para_as_value.c
  ${PROJECT_SOURCE_DIR}/src/math/asinh.c
  ${PROJECT_SOURCE_DIR}/src/math/e_acosh.c
  ${PROJECT_SOURCE_DIR}/src/math/e_atanh.c
  ${PROJECT_SOURCE_DIR}/src/math/factorial.c
  ${PROJECT_SOURCE_DIR}/src/math/fmax.c
  ${PROJECT_SOURCE_DIR}/src/math/fmin.c
  ${PROJECT_SOURCE_DIR}/src/math/isnan.c
  ${PROJECT_SOURCE_DIR}/src/math/s_log1p.c
  ${PROJECT_SOURCE_DIR}/src/solver/calc.c
  ${PROJECT_SOURCE_DIR}/src/solver/calc_event.c
  ${PROJECT_SOURCE_DIR}/src/solver/calc_initial_assignment.c
  ${PROJECT_SOURCE_DIR}/src/solver/calc_k.c
  ${PROJECT_SOURCE_DIR}/src/solver/calc_temp_value.c
  ${PROJECT_SOURCE_DIR}/src/solver/create_calc_object_list.c
  ${PROJECT_SOURCE_DIR}/src/solver/forwarding_value.c
  ${PROJECT_SOURCE_DIR}/src/solver/initialize_delay_val.c
  ${PROJECT_SOURCE_DIR}/src/solver/linear_approximation.c
  ${PROJECT_SOURCE_DIR}/src/solver/lu_decomposition.c
  ${PROJECT_SOURCE_DIR}/src/solver/lu_solve.c
  ${PROJECT_SOURCE_DIR}/src/solver/simulate_explicit.c
  ${PROJECT_SOURCE_DIR}/src/solver/simulate_implicit.c
  ${PROJECT_SOURCE_DIR}/src/solver/substitute_delay_val.c
  ${PROJECT_SOURCE_DIR}/src/util/chomp.c
  ${PROJECT_SOURCE_DIR}/src/util/dbg_printf.c
  ${PROJECT_SOURCE_DIR}/src/util/dupstr.c
  ${PROJECT_SOURCE_DIR}/src/util/get_end_cycle.c
  ${PROJECT_SOURCE_DIR}/src/util/my_time.c
  ${PROJECT_SOURCE_DIR}/src/util/prg_printf.c
  ${PROJECT_SOURCE_DIR}/src/util/set_seed.c
  ${PROJECT_SOURCE_DIR}/src/util/show_para.c
  ${PROJECT_SOURCE_DIR}/src/util/show_sp.c
  ${PROJECT_SOURCE_DIR}/src/util/str_is_number.c
  )

set (SOURCES_LIB ${SOURCES_COMMON} ${PROJECT_SOURCE_DIR}/src/lib_main.c)
set (SOURCES_MAIN ${SOURCES_COMMON} ${PROJECT_SOURCE_DIR}/src/main.c)
set (SOURCES_TEST test.c)

set (INSTALL_DOCS)
set (INSTALL_DOCS ${INSTALL_DOCS}
  ${PROJECT_SOURCE_DIR}/API.txt
  ${PROJECT_SOURCE_DIR}/AUTHORS.txt
  ${PROJECT_SOURCE_DIR}/COPYING.txt
  ${PROJECT_SOURCE_DIR}/FUNDING.txt
  ${PROJECT_SOURCE_DIR}/LICENSE.txt
  ${PROJECT_SOURCE_DIR}/LICENSE-dSFMT.txt
  ${PROJECT_SOURCE_DIR}/NEWS.txt
  ${PROJECT_SOURCE_DIR}/README.txt
  ${PROJECT_SOURCE_DIR}/VERSION.txt
  )

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
INSTALL(FILES ${INSTALL_DOCS} DESTINATION share/libsbmlsim)
# Install sample SBML, result files
FILE(GLOB sample_files "${PROJECT_SOURCE_DIR}/examples/sample.*")
INSTALL(FILES ${sample_files} DESTINATION share/libsbmlsim)
# Install C and C++ code
INSTALL(FILES ${PROJECT_SOURCE_DIR}/examples/c/test.c ${PROJECT_SOURCE_DIR}/examples/c/Makefile DESTINATION share/libsbmlsim/c)
INSTALL(FILES ${PROJECT_SOURCE_DIR}/examples/cpp/test.cpp ${PROJECT_SOURCE_DIR}/examples/cpp/Makefile DESTINATION share/libsbmlsim/cpp)
# Install dependent libraries on Windows
if(MSVC)
  FILE(GLOB windows_depend_dlls "${PROJECT_SOURCE_DIR}/depend/bin/*.dll")
  INSTALL(FILES ${windows_depend_dlls} DESTINATION bin)
  INSTALL(FILES ${PROJECT_SOURCE_DIR}/README.txt DESTINATION .)
endif()

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

if(WITH_CSHARP)
add_subdirectory(bindings/csharp)
endif()

# for installer
INCLUDE(InstallRequiredSystemLibraries)                                                                                                                                           

SET(CPACK_PACKAGE_DESCRIPTION "LibSBMLSim")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "LibSBMLSim: The library for simulating SBML models")
SET(CPACK_PACKAGE_NAME "libsbmlsim")
SET(CPACK_PACKAGE_VENDOR "LibSBMLSim development team")
SET(CPACK_PACKAGE_CONTACT "LibSBMLSim development Team <sbmlsim@fun.bio.keio.ac.jp>")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/README.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE.txt")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "1")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.txt")
SET(CPACK_RESOURCE_FILE_WELCOME "${PROJECT_SOURCE_DIR}/templates/Welcome.txt")

# for Source code archive
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")
SET(CPACK_SOURCE_GENERATOR "TGZ;ZIP") 
#SET(CPACK_INCLUDE_TOPLEVEL_DIRECTORY "1")
#SET(CPACK_SOURCE_TOPLEVEL_TAG "Source") 
#SET(CPACK_SOURCE_IGNORE_FILES "/CVS/;\\\\.DS_Store;/\\\\.svn/;/\\\\.git/;\\\\.swp$;\\\\.#;.*~;/templates/;/diff/;/simulation_results/;/bindings/;^l?v*;build_.*\\\\.sh;depend-.*\\\\.zip;ev_func_piece\\\\.xml") 
SET(CPACK_SOURCE_IGNORE_FILES 
  "/CVS/;/\\\\.svn/;/\\\\.git/;\\\\.DS_Store;\\\\.swp$"
  "/build/;^bindings/java;^bindings/perl;^bindings/python;/diff/;/simulation_results/"
  "testcases/simulateSBML"
  "testcases/cases/"
  "testcases/.*-results\\\\.csv"
  "testcases/out\\\\.csv"
  "build_csharp\\\\.sh;build_java\\\\.sh;build_python\\\\.sh;depend-.*\\\\.zip;result\\\\.zip"
  "Makefile\\\\.dist;ev_func_piece\\\\.xml;00remove\\\\.sh"
)

# 32bit or 64bit
IF(CMAKE_SIZEOF_VOID_P EQUAL 4)
  SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386")
  SET(CPACK_RPM_PACKAGE_ARCHITECTURE "i386")
  SET(PACKAGE_ARCH_EXT "x86")
ELSE()
  SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
  SET(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
  SET(PACKAGE_ARCH_EXT "x64")
ENDIF()

IF(WIN32 AND NOT UNIX)
  SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  #SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/Utilities/Release\\\\InstallIcon.bmp")
  #SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\MyExecutable.exe")
  SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}")
  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\fun.bio.keio.ac.jp\\\\software\\\\libsbmlsim")
  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\bio.keio.ac.jp\\\\software\\\\libsbmlsim")
  SET(CPACK_NSIS_CONTACT "sbmlsim@fun.bio.keio.ac.jp")
  SET(CPACK_NSIS_MODIFY_PATH ON)
ELSE(WIN32 AND NOT UNIX)
  #  SET(CPACK_INCLUDE_TOPLEVEL_DIRECTORY 0)
  SET(CPACK_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
  SET(CPACK_SET_DESTDIR "ON")
  SET(CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_VERSION}-${PACKAGE_ARCH_EXT})
  IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(CPACK_GENERATOR "DEB")
    SET(CPACK_INSTALL_PREFIX "/usr")
    SET(CPACK_RPM_PACKAGE_LICENSE "LGPL")
    SET(CPACK_RPM_PACKAGE_GROUP "Libraries/Development")
    SET(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
    SET(CPACK_DEBIAN_PACKAGE_SECTION "Libraries")
    STRING(TOLOWER "${CPACK_PACKAGE_NAME}" CPACK_PACKAGE_NAME_LOWERCASE) 
    SET(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME_LOWERCASE}_${CPACK_PACKAGE_VERSION}_${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")   
  ELSEIF(APPLE)
    SET (CPACK_GENERATOR "PackageMaker")
    IF(${CMAKE_SYSTEM_VERSION} MATCHES "^5\\.5") 
      SET(MACOSX_NAME "puma")
    ELSEIF(${CMAKE_SYSTEM_VERSION} MATCHES "^6\\.8") 
      SET(MACOSX_NAME "jaguar")
    ELSEIF(${CMAKE_SYSTEM_VERSION} MATCHES "^7\\.9") 
      SET(MACOSX_NAME "panther")
    ELSEIF(${CMAKE_SYSTEM_VERSION} MATCHES "^8\\.11") 
      SET(MACOSX_NAME "tiger")
    ELSEIF(${CMAKE_SYSTEM_VERSION} MATCHES "^9\\.8") 
      SET(MACOSX_NAME "leopard")
    ELSEIF(${CMAKE_SYSTEM_VERSION} MATCHES "^10\\.8") 
      SET(MACOSX_NAME "snowleopard")
    ELSEIF(${CMAKE_SYSTEM_VERSION} MATCHES "^11\\.4") 
      SET(MACOSX_NAME "lion")
    ELSEIF(${CMAKE_SYSTEM_VERSION} MATCHES "^12\\.0") 
      SET(MACOSX_NAME "mountainlion")
    ELSE()
      SET(MACOSX_NAME "unknown")
    ENDIF()
    SET(CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-macosx-${MACOSX_NAME}-${PACKAGE_ARCH_EXT})
  ENDIF()
  #SET(CPACK_STRIP_FILES "bin/MyExecutable")
  #SET(CPACK_SOURCE_STRIP_FILES "")
ENDIF(WIN32 AND NOT UNIX)

#SET(CPACK_PACKAGE_EXECUTABLES "MyExecutable" "My Executable")
INCLUDE(CPack)
