cmake_minimum_required(VERSION 2.6)
project(libsbmlsim)

# cmake variables
set(CMAKE_C_FLAGS "-g -O2 -Wall -ansi -pedantic")
# language bindings
option(WITH_JAVA     "Generate Java language bindings."   OFF)
option(WITH_PYTHON   "Generate Python language bindings." OFF)
option(WITH_PERL     "Generate Perl language bindings."   OFF)

# - Try to find libSBML
# Once done this will define
#  LIBSBML_FOUND - System has libSBML
#  LIBSBML_INCLUDE_DIRS - The libSBML include directories
#  LIBSBML_LIBRARIES - The libraries needed to use libSBML
#  LIBSBML_DEFINITIONS - Compiler switches required for using libSBML

find_package(PkgConfig)
pkg_check_modules(PC_LIBSBML QUIET libsbml)
set(LIBSBML_DEFINITIONS ${PC_LIBSBML_CFLAGS_OTHER})

find_path(LIBSBML_INCLUDE_DIR sbml/SBMLTypes.h
  HINTS ${PC_LIBSBML_INCLUDEDIR} ${PC_LIBSBML_INCLUDE_DIRS}
  PATH_SUFFIXES sbml )

find_library(LIBSBML_LIBRARY NAMES libsbml.lib sbml
  HINTS ${PC_LIBSBML_LIBDIR} ${PC_LIBSBML_LIBRARY_DIRS} )

set(LIBSBML_LIBRARIES ${LIBSBML_LIBRARY} )
set(LIBSBML_INCLUDE_DIRS ${LIBSBML_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBSBML_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LIBSBML DEFAULT_MSG LIBSBML_LIBRARY LIBSBML_INCLUDE_DIR)

add_subdirectory(src)

message(STATUS "
----------------------------------------------------------------------
libSBML version ${PACKAGE_VERSION}
----------------------------------------------------------------------

More information and the latest version are available online at
http://sbml.org/Software/libSBML

Please report problems using the issue tracker at
http://sbml.org/Software/libsbml/issue-tracker

To contact the developers directly, email libsbml-team@caltech.edu

Configured on host "$ENV{COMPUTERNAME}" ${HOSTNAME}
  host type                     = ${CMAKE_SYSTEM_NAME}
  host operating system         = ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}
  host CPU architecture         = ${CMAKE_SYSTEM_PROCESSOR}

General build flags:
  CC                            = ${CMAKE_C_COMPILER}
  CXX                           = ${CMAKE_CXX_COMPILER}
  CPPFLAGS                      = ${BUILD_DEFINITIONS}
  CFLAGS                        = ${CMAKE_C_FLAGS}
  CXXFLAGS                      = ${CMAKE_CXX_FLAGS}
  LDFLAGS                       = ${CMAKE_EXE_LINKER_FLAGS}

SBML parser library configuration:
  SBML library                   = ${LIBSBML_LIBRARY}
  CPPFLAGS for SBML              = -I${LIBSBML_INCLUDE_DIR}
  LIBS for SBML                  = ${LIBSBML_LIBRARY}

Other libSBMLSim configuration settings:
  Installation $prefix          = ${CMAKE_INSTALL_PREFIX}")

