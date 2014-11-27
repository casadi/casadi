cmake_minimum_required(VERSION 2.8.6)

find_package(SWIG REQUIRED)

# NOTE: you are not allowd to generate header files with configure_file in this module.
# See https://github.com/casadi/casadi/issues/646

# cmake_policy propagation rules
if(POLICY CMP0011)
cmake_policy(SET CMP0011 OLD)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)

# Include the use SWIG file
include(${CMAKE_MODULE_PATH}/UseSWIGMOD.cmake)

# Set misc SWIG flags
set(USE_SWIG_FLAGS ${USE_SWIG_FLAGS} "-DSWIG_TYPE_TABLE=casadi")
set(USE_SWIG_FLAGS ${USE_SWIG_FLAGS} "-fvirtual")
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../doc/doc.i)
  set(USE_SWIG_FLAGS ${USE_SWIG_FLAGS} "-DWITH_DOXDOC")
endif()

if(WITH_WERROR)
  # This will make the swig command return an error so that "make" will fail,
  # but swig will emit the output files anyways so subsequent calls to "make"
  # will succeed. At least buildbot will report failures.
  set(USE_SWIG_FLAGS ${USE_SWIG_FLAGS} "-Werror")
  set(MAYBE_WERROR "-Werror")
else()
  set(MAYBE_WERROR "")
endif()

if(WITH_PYTHON)
   add_subdirectory(python)
endif()

if(WITH_MATLAB)
   add_subdirectory(matlab)
endif()

if(WITH_JSON)
   add_subdirectory(json)
endif()
