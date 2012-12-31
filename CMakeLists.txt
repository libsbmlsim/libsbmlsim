cmake_minimum_required(VERSION 2.6)
project(libsbmlsim)

# cmake variables
if(MSVC)
  add_definitions(/W3)
else(MSVC)
  #  add_definitions(-Wall -Wextra -ansi -pedantic -fno-strict-aliasing)
  add_definitions(-Wall -Wextra -fno-strict-aliasing)
endif(MSVC)
set(CMAKE_C_FLAGS "-O2")

# build mode
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# language bindings
option(WITH_JAVA     "Generate Java language bindings."   OFF)
option(WITH_PYTHON   "Generate Python language bindings." OFF)
option(WITH_CSHARP   "Generate C# language bindings." OFF)
if(NOT WIN32)
  #  option(WITH_PERL   "Generate Perl language bindings."   OFF)
  option(WITH_RUBY   "Generate Ruby language bindings."   OFF)
endif(NOT WIN32)

# - Try to find libSBML
# Once done this will define
#  LIBSBML_FOUND - System has libSBML
#  LIBSBML_INCLUDE_DIRS - The libSBML include directories
#  LIBSBML_LIBRARIES - The libraries needed to use libSBML
#  LIBSBML_DEFINITIONS - Compiler switches required for using libSBML

find_package(PkgConfig)
pkg_check_modules(PC_LIBSBML QUIET libsbml)
set(LIBSBML_DEFINITIONS ${PC_LIBSBML_CFLAGS_OTHER})

find_path(LIBSBML_INCLUDE_DIR
  NAMES sbml/SBMLTypes.h
  PATHS /usr/include /usr/local/include /opt/local/include
        ${CMAKE_SOURCE_DIR}/depend/include
  HINTS ${PC_LIBSBML_INCLUDEDIR} ${PC_LIBSBML_INCLUDE_DIRS} )

find_library(LIBSBML_LIBRARY
  NAMES libsbml.lib sbml
  PATHS /usr/lib /usr/local/lib /opt/local/lib
        ${CMAKE_SOURCE_DIR}/depend/lib
  HINTS ${PC_LIBSBML_LIBDIR} ${PC_LIBSBML_LIBRARY_DIRS} )

find_library(LIBBZIP2_LIBRARY
  NAMES bzip2.lib bz2 libbz2.lib
  PATHS /usr/lib /usr/local/lib /opt/local/lib
        ${CMAKE_SOURCE_DIR}/depend/lib
  HINTS ${PC_LIBBZIP2} ${PC_LIBBZIP2_LIBRARY_DIRS} )

set(LIBSBML_LIBRARIES ${LIBSBML_LIBRARY} ${LIBBZIP2_LIBRARY} )
set(LIBSBML_INCLUDE_DIRS ${LIBSBML_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBSBML_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LIBSBML DEFAULT_MSG LIBSBML_LIBRARY LIBSBML_INCLUDE_DIR)

add_subdirectory(src)

message(STATUS "
----------------------------------------------------------------------
libSBMLSim version ${PACKAGE_VERSION}
----------------------------------------------------------------------

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

