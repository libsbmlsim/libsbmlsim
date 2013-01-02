# This is a CMake example for CSharp

cmake_minimum_required(VERSION 2.8)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

FIND_PACKAGE(SWIG 2.0 REQUIRED)
INCLUDE(${SWIG_USE_FILE})

# Only 32bit build on Windows
#if (WIN32)
#  set(CSHARP_PLATFORM "x64" CACHE STRING "C# target platfomr: x86, x64, anycpu, or itanium")
#endif (WIN32)

# Find C#
FIND_PACKAGE(CSharp REQUIRED)
INCLUDE(${CSHARP_USE_FILE})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src/libsbmlsim ${LIBSBML_INCLUDE_DIR})

# Make sure the nested directory structure exists
set(CSHARP_SOURCE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CSharpSources CACHE INTERNAL "")
set(CSHARP_BINARY_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CSharpBinaries CACHE INTERNAL "")
file(MAKE_DIRECTORY ${CSHARP_SOURCE_DIRECTORY})
file(MAKE_DIRECTORY ${CSHARP_BINARY_DIRECTORY})

# Create swig target
set(CMAKE_SWIG_OUTDIR ${CSHARP_SOURCE_DIRECTORY})
#SET(CMAKE_SWIG_FLAGS -namespace \"libsbmlsim\" ${CMAKE_SWIG_GLOBAL_FLAGS} ${CMAKE_SWIG_FLAGS})
SET(CMAKE_SWIG_FLAGS -namespace \"libsbmlsimcs\" ${CMAKE_SWIG_GLOBAL_FLAGS} ${CMAKE_SWIG_FLAGS})
#if (UNIX)
#  set(CMAKE_SWIG_FLAGS -dllimport \"libsbmlsim\")
#else (WIN32)
#  set(CMAKE_SWIG_FLAGS -dllimport \"sbmlsim\")
#endif (UNIX)
set(CMAKE_SWIG_FLAGS -dllimport \"libsbmlsimcs\")

SET_SOURCE_FILES_PROPERTIES(libsbmlsim.i PROPERTIES CPLUSPLUS OFF)
SET_SOURCE_FILES_PROPERTIES(libsbmlsim.i PROPERTIES SWIG_FLAGS "-includeall")
SWIG_ADD_MODULE(libsbmlsimcs csharp libsbmlsim.i ${SOURCES_LIB})
SWIG_LINK_LIBRARIES(libsbmlsimcs ${LIBSBML_LIBRARIES})

set_target_properties(libsbmlsimcs PROPERTIES OUTPUT_NAME "libsbmlsimcs")
#set_target_properties(libsbmlsimcs PROPERTIES OUTPUT_NAME "libsbmlsim")
set_target_properties(libsbmlsimcs PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CSHARP_BINARY_DIRECTORY})
if ( UNIX )
  set_target_properties(libsbmlsimcs PROPERTIES PREFIX "")
  if (APPLE)
    set_target_properties (libsbmlsimcs PROPERTIES SUFFIX ".bundle")
  else (APPLE)
    set_target_properties(libsbmlsimcs PROPERTIES SUFFIX ".so")
  endif()
else ( WIN32 )
  set_target_properties(libsbmlsimcs PROPERTIES PREFIX "")
  set_target_properties(libsbmlsimcs PROPERTIES SUFFIX ".dll")
  foreach ( CMAKE_CONFIGURATION_TYPE ${CMAKE_CONFIGURATION_TYPES} )
    string(TOUPPER ${CMAKE_CONFIGURATION_TYPE} CMAKE_CONFIGURATION_TYPE)
    set_target_properties(libsbmlsimcs
      PROPERTIES LIBRARY_OUTPUT_DIRECTORY_${CMAKE_CONFIGURATION_TYPE} "${CSHARP_BINARY_DIRECTORY}")
    set_target_properties(libsbmlsimcs
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY_${CMAKE_CONFIGURATION_TYPE} "${CSHARP_BINARY_DIRECTORY}")
  endforeach( )
endif( UNIX )

# Configure AssemblyInfo.cs
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/AssemblyInfo.cs.in
  ${CSHARP_SOURCE_DIRECTORY}/AssemblyInfo.cs
  )

# Add managed wrapper
csharp_add_library(
  libsbmlsimcsP
  ${CSHARP_SOURCE_DIRECTORY}/*.cs
  )
add_dependencies(libsbmlsimcsP libsbmlsimcs)

file(COPY ${PROJECT_SOURCE_DIR}/examples/csharp/Test.cs DESTINATION ${CSHARP_BINARY_DIRECTORY})
file(COPY ${PROJECT_SOURCE_DIR}/examples/csharp/README.txt DESTINATION ${CSHARP_BINARY_DIRECTORY})
file(COPY ${PROJECT_SOURCE_DIR}/examples/sample.xml DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../)

# Install targets
INSTALL(FILES ${PROJECT_SOURCE_DIR}/examples/csharp/Test.cs ${PROJECT_SOURCE_DIR}/examples/csharp/README.txt ${CSHARP_BINARY_DIRECTORY}/libsbmlsimcsP.dll DESTINATION share/libsbmlsim/csharp)
INSTALL(TARGETS libsbmlsimcs 
  LIBRARY DESTINATION share/libsbmlsim/csharp
  ARCHIVE DESTINATION share/libsbmlsim/csharp
  )

