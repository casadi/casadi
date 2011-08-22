cmake_minimum_required(VERSION 2.0)
project(opticon C CXX)

# Options:
option(ENABLE_STATIC "Compile static libraries (.a on Linux)" ON)
option(ENABLE_SHARED "Compile shared libraries (.so on Linux)" OFF)
option(WITH_SUPERLU "Compile the SuperLU interface and source" ON)
option(WITH_MUSCOD "Compile the MUSCOD-II interface" OFF)
# option(WITH_IPOPT "Compile the Ipopt interface" ON)
# option(WITH_ACADO "Compile the ACADO interface" ON)
# option(WITH_SUNDIALS "Compile the Sundials interface" ON)
option(WITH_JIT "Enable just-in-time compilation" OFF)
option(WITH_DOC "Enable documentation generation" OFF)
option(WITH_PYTHON_INTERRUPTS "With interrupt handling inside python interface" OFF)
option(WITH_PRINTME "With support for printme operator (might slow things down)" OFF)
option(WITH_GSL "Compile the GSL interface" ON)
option(WITH_OPENMP "Compile with parallelization support" OFF)
option(WITH_OOQP "Enable OOQP interface" OFF)
option(WITH_FORTRAN "Enable Fortran linking, if this is set to OFF, no Fortran linking will be invoked (should not be needed, CMakeDetermineFortranCompiler should do the job)" ON)

if(WITH_FORTRAN)
  # Enable Fortran linking (needed e.g. for LAPACK/BLAS), if Fortran is available
  include(CMakeDetermineFortranCompiler)
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
  endif(CMAKE_Fortran_COMPILER)
endif(WITH_FORTRAN)

# For code optimization
if(CMAKE_BUILD_TYPE)
else(CMAKE_BUILD_TYPE)
SET(CMAKE_BUILD_TYPE Release)
# SET(CMAKE_BUILD_TYPE Debug)
endif(CMAKE_BUILD_TYPE)

# Necessary for Visual C++
add_definitions(-D_USE_MATH_DEFINES)
add_definitions(-D_SCL_SECURE_NO_WARNINGS)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

if(WITH_JIT)
add_definitions(-DWITH_JIT)
endif(WITH_JIT)

if(WITH_PRINTME)
add_definitions(-DWITH_PRINTME)
endif(WITH_PRINTME)

include_directories(.)

# set(CMAKE_VERBOSE_MAKEFILE 0)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi -pedantic -Wall -Wno-sign-compare")

# the following commands are needed to fix a problem with the libraries
# for linux 64 bits
if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
  message(STATUS "x86_64 architecture detected - setting flag -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
endif("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

# enabling openmp support if requesed
if(WITH_OPENMP)
find_package (OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif(OPENMP_FOUND)
endif(WITH_OPENMP)

#SET(LIBRARY_OUTPUT_PATH  ${CMAKE_CURRENT_SOURCE_DIR}/lib)
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/)

add_subdirectory(casadi)
add_subdirectory(optimal_control)
add_subdirectory(toolboxes)
add_subdirectory(examples)
add_subdirectory(external_packages)
add_subdirectory(interfaces)
add_subdirectory(experimental/greg)
add_subdirectory(experimental/joel)

if(WITH_MUSCOD)
add_subdirectory(muscod_interface)
endif(WITH_MUSCOD)

find_package (SWIG)
if(SWIG_FOUND)
add_subdirectory(swig EXCLUDE_FROM_ALL) # the swig interfaces are not compiled when typing "make"
add_custom_target(python DEPENDS _casadi) # the python interface is created with the command "make python"
add_custom_target(octave DEPENDS casadi_octave) # the octave interface is created with the command "make octave"
endif(SWIG_FOUND)

add_subdirectory(documentation/tutorials EXCLUDE_FROM_ALL)
add_custom_target(tutor)

if(WITH_DOC)
set(DOXYFILE_LATEX "NO")
set(DOXYFILE_IN_DIR ${CMAKE_SOURCE_DIR}/documentation/api-doc)
set(DOXYFILE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/casadi)
include(UseDoxygen OPTIONAL)
endif(WITH_DOC)

