cmake_minimum_required(VERSION 2.6)
project(casadi_matlab_interface C CXX)

include_directories(/usr/local/matlab/extern/include)

set(CASADI_MATLAB_SRCS
  matlab_interface.h
  matlab_interface.cxx
)
add_library(casadi_matlab SHARED ${CASADI_MATLAB_SRCS})

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../../cmake_modules)

find_package (SWIG REQUIRED QUIET)
include(${SWIG_USE_FILE})


FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Find SWIG and include the use swig file
# include (../../../cmake_modules/UseSWIGMOD.cmake)

# set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} ${CMAKE_DL_LIBS})

# set(CASADI_SWIG_MAIN_MODULE casadi_main)

set_source_files_properties (test_swig_matlab.i CPLUSPLUS ON)
set_source_files_properties (test_swig_matlab.i PROPERTIES SWIG_FLAGS "-includeall")
set(CMAKE_SWIG_FLAGS "")
# set(USE_SWIG_FLAGS ${USE_SWIG_FLAGS} "${CASADI_SWIG_FLAGS}")

# PYTHON interface
SWIG_ADD_MODULE (test_swig_matlab python test_swig_matlab.i test_swig_matlab.cpp)

# add_dependencies(_casadi casadi) #workaround to add dependencies of casadi to _casadi (SWIG_ADD_MODULE does not add dependencies from c++ source files)
# SWIG_LINK_LIBRARIES (casadi_matlab ${PYTHON_LIBRARIES}  ${SWIG_LINK_LIBRARIES_LIST} casadi_matlab )

