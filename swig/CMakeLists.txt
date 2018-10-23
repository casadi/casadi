cmake_minimum_required(VERSION 2.8.6)

if (NOT SWIG_IMPORT)
  find_package(SWIG REQUIRED)
endif()

if (WITH_EXTENDING_CASADI)
   add_subdirectory(extending_casadi)
endif()

# NOTE: you are not allowed to generate header files with configure_file in this module.
# See https://github.com/casadi/casadi/issues/646

# cmake_policy propagation rules
if(POLICY CMP0011)
cmake_policy(SET CMP0011 OLD)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)

if(WITH_WERROR)
  # This will make the swig command return an error so that "make" will fail,
  # but swig will emit the output files anyways so subsequent calls to "make"
  # will succeed. At least buildbot will report failures.
  set(USE_SWIG_FLAGS ${USE_SWIG_FLAGS} "-Werror")
  set(MAYBE_WERROR "-Werror")
else()
  set(MAYBE_WERROR "")
endif()

# Include the use SWIG file
#include(UseSWIG)
include(${CMAKE_MODULE_PATH}/UseSWIGMOD.cmake)

# Set SWIG flags
set(CMAKE_SWIG_FLAGS $ENV{EXTRA_SWIG_FLAGS})
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DSWIG_TYPE_TABLE=casadi")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-fvirtual")
#set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-fcompact")
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../doc/doc.i)
  set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DWITH_DOXDOC")
endif()

if(WITH_DEPRECATED_FEATURES)
  set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DWITH_DEPRECATED_FEATURES")
endif()

if(WITH_PYTHON)
   add_subdirectory(python)
endif()

if(WITH_MATLAB OR WITH_OCTAVE)
   add_subdirectory(matlab)
endif()

if(WITH_JSON)
   add_subdirectory(json)
endif()

install(FILES casadi.i DESTINATION ${INCLUDE_PREFIX}/casadi)
install(FILES doc.i DESTINATION ${INCLUDE_PREFIX}/casadi)
install(FILES doc_merged.i DESTINATION ${INCLUDE_PREFIX}/casadi)
install(FILES python/casadi_numpy.hpp DESTINATION ${INCLUDE_PREFIX}/casadi)
install(FILES ../test/internal/valgrind-casadi.supp DESTINATION ${INCLUDE_PREFIX}/casadi)
install(FILES ../test/internal/valgrind-python.supp DESTINATION ${INCLUDE_PREFIX}/casadi)
