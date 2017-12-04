# This is a CMake example for Python

cmake_minimum_required(VERSION 2.8)

FIND_PACKAGE(SWIG 2.0 REQUIRED)
INCLUDE(${SWIG_USE_FILE})

INCLUDE(../../../cmake/FindPythonLibsNew.cmake)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src/libsbmlsim ${LIBSBML_INCLUDE_DIR})

SET(CMAKE_SWIG_FLAGS "")

SET_SOURCE_FILES_PROPERTIES(libsbmlsim.i PROPERTIES CPLUSPLUS OFF)
SET_SOURCE_FILES_PROPERTIES(libsbmlsim.i PROPERTIES SWIG_FLAGS "-includeall")
IF(${CMAKE_VERSION} VERSION_LESS "3.8.0")
  SWIG_ADD_MODULE(libsbmlsim python libsbmlsim.i ${SOURCES_LIB})
ELSE()
  SWIG_ADD_LIBRARY(libsbmlsim
    LANGUAGE python
    SOURCES libsbmlsim.i ${SOURCES_LIB})
ENDIF()
SWIG_LINK_LIBRARIES(libsbmlsim ${PYTHON_LIBRARIES} ${LIBSBML_LIBRARIES})

# show messages
MESSAGE( STATUS "PYTHON_PREFIX: " ${PYTHON_PREFIX} )               #  - path to the Python installation
MESSAGE( STATUS "PYTHON_LIBRARIES: " ${PYTHON_LIBRARIES} )         #  - path to the python library
MESSAGE( STATUS "PYTHON_INCLUDE_DIRS: " ${PYTHON_INCLUDE_DIRS} )   #  - path to where Python.h is found
MESSAGE( STATUS "PYTHON_SITE_PACKAGES: " ${PYTHON_SITE_PACKAGES} ) #  - path to installation site-packages
MESSAGE( STATUS "PYTHON_IS_DEBUG: " ${PYTHON_IS_DEBUG} )           #  - whether the Python interpreter is a debug build

file(COPY ${PROJECT_SOURCE_DIR}/examples/python/Test.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/examples/python/README.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/examples/sample.xml DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../)

# Install targets
INSTALL(FILES ${PROJECT_SOURCE_DIR}/examples/python/Test.py
  ${PROJECT_SOURCE_DIR}/examples/python/gridSearch.py
  ${PROJECT_SOURCE_DIR}/examples/python/data.csv
  ${PROJECT_SOURCE_DIR}/examples/python/simple.xml
  ${PROJECT_SOURCE_DIR}/examples/python/README.txt DESTINATION share/libsbmlsim/python)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/libsbmlsim.py DESTINATION ${PYTHON_SITE_PACKAGES})

INSTALL(TARGETS _libsbmlsim
  LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES}
  ARCHIVE DESTINATION ${PYTHON_SITE_PACKAGES}
  )

