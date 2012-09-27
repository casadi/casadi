cmake_minimum_required(VERSION 2.0)
project(opticon C CXX)

# find the macports python libs if they exist
if(APPLE)
  set(CMAKE_LIBRARY_PATH /opt/local/lib ${CMAKE_LIBRARY_PATH})
endif()

# What is this, why is is needed and why does it not appear to work?
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif()

# We use the FeatureSummary module to get a more readable CMake output
include(FeatureSummary)

# Options:
option(ENABLE_STATIC "Compile static libraries (.a on Linux)" ON)
option(ENABLE_SHARED "Compile shared libraries (.so on Linux)" OFF)
option(WITH_MUSCOD "Compile the MUSCOD-II interface" OFF)
option(WITH_DL "Enable dynamic loading of functions" OFF)
option(WITH_LLVM "Use LLVM for just-in-time compilation" OFF)
option(WITH_DOC "Enable documentation generation" OFF)
option(WITH_PYTHON_INTERRUPTS "With interrupt handling inside python interface" OFF)
# option(WITH_GSL "Compile the GSL interface" ON)
option(WITH_OPENMP "Compile with parallelization support" OFF)
option(WITH_OOQP "Enable OOQP interface" OFF)
option(WITH_FORTRAN "Enable Fortran linking, if this is set to OFF, no Fortran linking will be invoked (should not be needed, CMakeDetermineFortranCompiler should do the job)" ON)
option(WITH_SWIG_SPLIT "Split SWIG wrapper generation into multiple modules" OFF) 
option(WITH_WORHP "Compile the WORHP interface" OFF) 
option(WITH_ACADO "Compile the interfaces to ACADO and qpOASES, if it can be found" OFF) 
option(WITH_EIGEN3 "Use Eigen3 backend for matrix operations" OFF) 
option(WITH_SUNDIALS "Compile the interface to Sundials (the source code for Sundials 2.5 is included)" ON)
option(WITH_QPOASES "Compile the interface to qpOASES (the source code for qpOASES 3.0beta is included)" ON)
option(WITH_CSPARSE "Compile the interface to CSparse (the source code for CSparse is included)" ON)

# Enable Fortran linking (needed e.g. for LAPACK/BLAS), if Fortran is available
if(WITH_FORTRAN)
  include(CMakeDetermineFortranCompiler)
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
  endif()
endif()

# For code optimization
if(CMAKE_BUILD_TYPE)
else(CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  # set(CMAKE_BUILD_TYPE Debug)
endif(CMAKE_BUILD_TYPE)

# Necessary for Visual C++
add_definitions(-D_USE_MATH_DEFINES)
add_definitions(-D_SCL_SECURE_NO_WARNINGS)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

if(WITH_DL)
  add_definitions(-DWITH_DL)
endif()

if(WITH_PRINTME)
  add_definitions(-DWITH_PRINTME)
endif()

include_directories(.)

# Check if the compiler supports C++11
foreach(CXX11FLAG_TEST "" "-std=c++11" "-std=c++0x" "-std=gnu++0x")
  if(NOT CXX11FLAG)
    # Try to compile test with the proposed flag
    try_compile(CXX11_SUCCESS
      ${CMAKE_BINARY_DIR} 
      ${CMAKE_SOURCE_DIR}/cmake_modules/check_cxx11.cpp
      CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${CXX11FLAG_TEST}
      OUTPUT_VARIABLE OUTPUT)
    
    # Save flag, if successful
    if(CXX11_SUCCESS)
      set(CXX11FLAG ${CXX11FLAG_TEST})
    endif()
  endif()
endforeach()

# If found, set flag
if(CXX11FLAG)
  message(STATUS "Flag for enabling C++11 features: ${CXX11FLAG}")
  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${CXX11FLAG})
  add_definitions(-DUSE_CXX11)
endif()

add_feature_info(using-c++11 USE_CXX11 "Using C++11 features (improves efficiency and is required for some examples).")

# set(CMAKE_VERBOSE_MAKEFILE 0)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi -pedantic -Wall -Wno-sign-compare")

# the following commands are needed to fix a problem with the libraries
# for linux 64 bits
if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
  message(STATUS "x86_64 architecture detected - setting flag -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
endif("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

# Use the included Sundials
if(WITH_SUNDIALS)
  set(SUNDIALS_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external_packages/sundials-2.5mod)
  set(SUNDIALS_LIBRARIES sundials-2.5mod)
endif()
add_feature_info(sundials-interface WITH_SUNDIALS "Interface to the ODE/DAE integrator suite SUNDIALS.")

# Use the included qpOASES
if(WITH_QPOASES)
  set(QPOASES_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external_packages/qpOASES/include)
  set(QPOASES_LIBRARIES qpoases)
  set(QPOASES_FOUND true)
endif()
add_feature_info(qpoases-interface WITH_QPOASES "Interface to the active-set QP solver qpOASES.")

# Use the included qpOASES
if(WITH_CSPARSE)
  set(CSPARSE_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external_packages/CSparse/Include)
  set(CSPARSE_LIBRARIES csparse)
endif()
add_feature_info(csparse-interface WITH_CSPARSE "Interface to the sparse direct linear solver CSparse.")

# enabling openmp support if requesed
if(WITH_OPENMP)
  find_package (OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DWITH_OPENMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -DWITH_OPENMP")
  endif()
endif()

#SET(LIBRARY_OUTPUT_PATH  ${CMAKE_CURRENT_SOURCE_DIR}/lib)
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/)

# Optional core dependencies
if(WITH_LLVM)
  # Core depends on LLVM for Just-in-time compilation
  find_package(LLVM REQUIRED)
  find_package(Threads REQUIRED)
  set(CASADI_DEPENDENCIES ${LLVM_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CMAKE_DL_LIBS})
  add_definitions(-DWITH_LLVM)
  include_directories( ${LLVM_INCLUDE_DIR} )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LLVM_CFLAGS}") # Note LLVM_CFLAGS, not LLVM_CXXFLAGS due to -fno_exceptions
else(WITH_LLVM)
  if(WITH_DL)
    # Core needs support for dynamic linking
    set(CASADI_DEPENDENCIES ${CMAKE_DL_LIBS})
  else(WITH_DL)
    # No dependencies
    set(CASADI_DEPENDENCIES)
  endif(WITH_DL)
endif(WITH_LLVM)
add_feature_info(just-in-time WITH_LLVM "Just-in-time compiliation via the LLVM compiler framework.")

# Optional auxillary dependencies
find_package(BLAS QUIET)

find_package(LAPACK)

find_package(IPOPT)
add_feature_info(ipopt-interface IPOPT_FOUND "Interface to the NLP solver Ipopt.")

if(WITH_ACADO)
  find_package(ACADO)
endif()

find_package(KNITRO)
add_feature_info(knitro-interface KNITRO_FOUND "Interface to the NLP solver KNITRO.")

find_package(LIFTOPT)
add_feature_info(liftopt-interface LIFTOPT_FOUND "Interface to the Lifted Newton implementation LiftOpt.")

find_package(CPLEX)
add_feature_info(cplex-interface CPLEX_FOUND "Interface to the QP solver CPLEX.")

# if(WITH_GSL)
#   find_package (GSL)
# endif()

if(WITH_OOQP)
  find_package (OOQP)
endif()
add_feature_info(ooqp-interface OOQP_FOUND "Interface to the QP solver OOQP.")

if(WITH_WORHP)
  find_package (WORHP)
endif()
add_feature_info(worhp-inteface WORHP_FOUND "Interface to the NLP solver Worhp.")

add_subdirectory(casadi)
add_subdirectory(optimal_control)
add_subdirectory(nonlinear_programming)
add_subdirectory(integration)
add_subdirectory(toolboxes)
add_subdirectory(examples)
add_subdirectory(external_packages)
add_subdirectory(interfaces)
add_subdirectory(experimental/greg EXCLUDE_FROM_ALL)
add_subdirectory(experimental/joel EXCLUDE_FROM_ALL)
#add_subdirectory(experimental/andrew EXCLUDE_FROM_ALL)

if(WITH_MUSCOD)
  add_subdirectory(muscod_interface)
endif(WITH_MUSCOD)

find_package(SWIG)
  if(SWIG_FOUND)
    find_package(MPI QUIET)
    add_subdirectory(swig EXCLUDE_FROM_ALL) # the swig interfaces are not compiled when typing "make"
    if(WITH_SWIG_SPLIT)
      add_custom_target(python DEPENDS _casadi_main _casadi_primitive _casadi_primitive_tools _casadi_noncore) # the python interface is created with the command "make python"
      add_custom_target(octave DEPENDS casadi_main_octave casadi_primitive_octave casadi_primitive_tools_octave casadi_noncore_octave) # the octave interface is created with the command "make octave"
    else(WITH_SWIG_SPLIT)
      add_custom_target(python DEPENDS _casadi) # the python interface is created with the command "make python"
      add_custom_target(octave DEPENDS casadi_octave ) # the octave interface is created with the command "make octave"
  endif(WITH_SWIG_SPLIT)
endif(SWIG_FOUND)

add_subdirectory(documentation/tutorials EXCLUDE_FROM_ALL)
add_custom_target(tutor)

if(WITH_DOC)
  set(DOXYFILE_LATEX "NO")
  set(DOXYFILE_IN_DIR ${CMAKE_SOURCE_DIR}/documentation/api-doc)
  set(DOXYFILE_OUTPUT_DIR  ${CMAKE_SOURCE_DIR}/documentation/api-doc)
  set(DOXYFILE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/casadi)
  include(UseDoxygen OPTIONAL)
endif()

add_subdirectory(documentation/examples/ctemplate)

# Print a summary
feature_summary(WHAT ALL)

