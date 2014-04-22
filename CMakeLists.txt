cmake_minimum_required(VERSION 2.8.6)
project(casadi C CXX)
set(PACKAGE_NAME "CasADi")
set(PACKAGE_VERSION "1.9.0+")
set(PACKAGE_DESCRIPTION "CasADi is a symbolic framework for automatic differentation and numeric optimization")
set(PACKAGE_AUTHOR "Joel Andersson")
set(PACKAGE_MAINTAINER "Joris Gillis")
set(PACKAGE_URL "http://casadi.org")

# For systems with git installed, find out revsion and description.
execute_process(COMMAND git rev-parse HEAD
OUTPUT_VARIABLE git_revision
OUTPUT_STRIP_TRAILING_WHITESPACE
ERROR_QUIET
)

execute_process(COMMAND git describe HEAD
OUTPUT_VARIABLE git_describe
OUTPUT_STRIP_TRAILING_WHITESPACE
ERROR_QUIET
)

# special version number for pkg-config, use n.n.n.n where n are all numbers
string(REGEX REPLACE "-g.*" "" PACKAGE_VERSION_PKGCONFIG "${git_describe}")
string(REPLACE "+-" "." PACKAGE_VERSION_PKGCONFIG "${PACKAGE_VERSION_PKGCONFIG}")

string(REPLACE "+-" "+" git_describe "${git_describe}")
string(REPLACE "-g" "." git_describe "${git_describe}")

set(PACKAGE_VERSION_FULL "${git_describe}")

# overwrite version numbers for official release
if((PACKAGE_VERSION_FULL STREQUAL "") OR (NOT PACKAGE_VERSION MATCHES "\\+"))
  set(PACKAGE_VERSION_PKGCONFIG "${PACKAGE_VERSION}")
  set(PACKAGE_VERSION_FULL "${PACKAGE_VERSION}")
endif()


# begin workaround for
# http://public.kitware.com/Bug/print_bug_page.php?bug_id=9220
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(language_support)
workaround_9220(Fortran Fortran_language_works)
if(Fortran_language_works)
  enable_language(Fortran OPTIONAL)
endif()
#end workaround

include(CppLint)
include(SpellCheck)
include(embed_resources)
include (GenerateExportHeader)

# find the macports python libs if they exist
if(APPLE)
  set(CMAKE_LIBRARY_PATH /opt/local/lib ${CMAKE_LIBRARY_PATH})
endif()

# Enforce cmake2.4+ style behaviour for
# mixed filename/full path linker directives
if(COMMAND cmake_policy)
cmake_policy(SET CMP0003 NEW)
endif()

# cmake_policy propagation rules
if(COMMAND cmake_policy)
if(POLICY CMP0011)
cmake_policy(SET CMP0011 OLD)
endif()
endif()

# cmake_policy avoid escaping in add_definitions
if(COMMAND cmake_policy)
cmake_policy(SET CMP0005 OLD)
endif()

# We use the FeatureSummary module to get a more readable CMake output
include(FeatureSummary)
if(NOT COMMAND add_feature_info)
  function(add_feature_info _name _enabled _description)
    message(STATUS "Detected: ${_name} ${_description}")
  endfunction()
  function(feature_summary)
  endfunction()
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  option(ENABLE_STATIC "Compile static libraries (.a on Linux)" ON)
  option(ENABLE_SHARED "Compile shared libraries (.so on Linux)" OFF)
else()
  option(ENABLE_STATIC "Compile static libraries (.a on Linux)" OFF)
  option(ENABLE_SHARED "Compile shared libraries (.so on Linux)" ON)
endif()
  
# Options:
if(WIN32)
  option(WITH_DL "Enable dynamic loading of functions" OFF)
else()
  option(WITH_DL "Enable dynamic loading of functions" ON)
endif()
option(WITH_DOC "Enable documentation generation" OFF)
option(ENABLE_EXPORT_ALL "Export all symbols to a shared library" OFF)

option(WITH_PYTHON_INTERRUPTS "With interrupt handling inside python interface" OFF)
# option(WITH_GSL "Compile the GSL interface" ON)
option(WITH_OPENMP "Compile with parallelization support" OFF)
option(WITH_OOQP "Enable OOQP interface" ON)
option(WITH_SQIC "Enable SQIC interface" OFF)
option(WITH_SLICOT "Enable SLICOT interface" OFF)
option(WITH_WORHP "Compile the WORHP interface" OFF)
option(WITH_SNOPT "Compile the SNOPT interface" OFF)
option(WITH_ACADO "Compile the interfaces to ACADO and qpOASES, if it can be found" OFF)
option(WITH_EIGEN3 "Use Eigen3 backend for matrix operations" OFF)
option(WITH_SUNDIALS "Compile the interface to Sundials (the source code for Sundials 2.5 is included)" ON)
option(WITH_QPOASES "Compile the interface to qpOASES (the source code for qpOASES 3.0beta is included)" ON)
option(WITH_DSDP "Compile the interface to DSDP (the source code for DSDP is included)" ON)
option(WITH_CSPARSE "Compile the interface to CSparse (the source code for CSparse is included)" ON)
option(WITH_IPOPT "Compile the interface to IPOPT" ON)
option(WITH_KNITRO "Compile the interface to KNITRO" ON)
option(WITH_CPLEX "Compile the interface to CPLEX" ON)
option(WITH_LAPACK "Compile the interface to LAPACK" ON)
option(WITH_OPENCL "Compile with OpenCL support" OFF)
option(WITH_PROFILING "Enable a built-in profiler to be switched used" OFF)
option(WITH_COVERAGE "Create coverage report" OFF)

# For code optimization
if(CMAKE_BUILD_TYPE)
else()
  set(CMAKE_BUILD_TYPE Release)
  #set(CMAKE_BUILD_TYPE Debug)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} UPPER_CMAKE_BUILD_TYPE)

# Necessary for Visual C++
add_definitions(-D_USE_MATH_DEFINES)
add_definitions(-D_SCL_SECURE_NO_WARNINGS)

if(WITH_DL)
  add_definitions(-DWITH_DL)
endif()
add_feature_info(dynamic-loading WITH_DL "Compile with support for dynamic loading of generated functions (needed for ExternalFunction)")

if(WITH_PRINTME)
  add_definitions(-DWITH_PRINTME)
endif()

include_directories(.)
include_directories(${CMAKE_BINARY_DIR})

# Check if the compiler supports C++11
foreach(CXX11FLAG_TEST "" "-std=gnu++11" "-std=c++11" "-std=gnu++0x" "-std=c++0x" "-std=c++11 --DUSE_TR1_HASHMAP" "-std=c++0x -DUSE_TR1_HASHMAP")
  if(NOT USE_CXX11)
    # Try to compile test with the proposed flag
    try_compile(USE_CXX11
      ${CMAKE_BINARY_DIR}
      ${CMAKE_MODULE_PATH}/check_cxx11.cpp
      CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${CXX11FLAG_TEST}
      OUTPUT_VARIABLE OUTPUT)

    # Save flag, if successful
    if(USE_CXX11)
      set(CXX11FLAG ${CXX11FLAG_TEST} CACHE INTERNAL "C++11 compile flag")
    endif()
  endif()
endforeach()

# If found, set flag
if(USE_CXX11)
  if(CXX11FLAG)
    message(STATUS "Flag needed for enabling C++11 features: ${CXX11FLAG}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX11FLAG}")
  else()
    message(STATUS "No flag needed for enabling C++11 features.")
  endif()
  add_definitions(-DUSE_CXX11)
endif()
add_feature_info(using-c++11 USE_CXX11 "Using C++11 features (improves efficiency and is required for some examples).")

if(CXX11FLAG)
  try_compile(HAS_COPYSIGN
    ${CMAKE_BINARY_DIR}
    ${CMAKE_MODULE_PATH}/check_copysign.cpp
    CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${CXX11FLAG}
    OUTPUT_VARIABLE OUTPUT)
  if(HAS_COPYSIGN)
    add_definitions(-DHAS_COPYSIGN)
  endif()

  try_compile(HAS_ERF
    ${CMAKE_BINARY_DIR}
    ${CMAKE_MODULE_PATH}/check_erf.cpp
    CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${CXX11FLAG}
    OUTPUT_VARIABLE OUTPUT)
  if(HAS_ERF)
    add_definitions(-DHAS_ERF)
  endif()

  try_compile(HAS_ERFINV
    ${CMAKE_BINARY_DIR}
    ${CMAKE_MODULE_PATH}/check_erfinv.cpp
    CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${CXX11FLAG}
    OUTPUT_VARIABLE OUTPUT)
  if(HAS_ERFINV)
    add_definitions(-DHAS_ERFINV)
  endif()
endif()


# set(CMAKE_VERBOSE_MAKEFILE 0)

# add compiler warnings
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
  "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR
  "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare -Wno-delete-non-virtual-dtor")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # I don't know the MSVC warnings, so do nothing for now
endif()

# the following commands are needed to fix a problem with the libraries
# for linux 64 bits
if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
  message(STATUS "x86_64 architecture detected - setting flag -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
endif()

if(WITH_COVERAGE)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fprofile-arcs -ftest-coverage")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage")
set(CMAKE_EXE_LINKER_FLAGS  " -lgcov -fprofile-arcs --coverage ${CMAKE_EXE_LINKER_FLAGS}")
endif()

if(MINGW)
  # Circumventing a bug in MinGW g++ v4.7.2, evoked by 752fa89355ffa
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-ipa-cp-clone")
endif()

# Use the included Sundials
if(WITH_SUNDIALS)
  set(SUNDIALS_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external_packages/sundials-2.5mod)
  set(SUNDIALS_LIBRARIES casadi_sundials)
endif()
add_feature_info(sundials-interface WITH_SUNDIALS "Interface to the ODE/DAE integrator suite SUNDIALS.")

# Use the included csparse
if(WITH_CSPARSE)
  set(CSPARSE_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external_packages/CSparse/Include)
  set(CSPARSE_LIBRARIES casadi_csparse)
endif()
add_feature_info(csparse-interface WITH_CSPARSE "Interface to the sparse direct linear solver CSparse.")

# TinyXML
set(TINYXML_LIBRARIES casadi_tinyxml)

# enabling openmp support if requesed
if(WITH_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DWITH_OPENMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -DWITH_OPENMP")
  endif()
endif()

#set(LIBRARY_OUTPUT_PATH  ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/)

# OpenCL
if(WITH_OPENCL)
  # Core depends on OpenCL for GPU calculations
  find_package(OpenCL REQUIRED)
  add_definitions(-DWITH_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIRS})
endif()
add_feature_info(opencl-support WITH_OPENCL "Enable just-in-time compiliation to CPUs and GPUs with OpenCL.")

# Optional auxillary dependencies
find_package(BLAS QUIET)
find_package(HSL QUIET)
find_package(WSMP QUIET)
find_package(METIS QUIET)
find_package(LibXml2 QUIET)
if(WITH_LAPACK)
  find_package(LAPACK)
endif()
add_feature_info(lapack-interface LAPACK_FOUND "Interface to LAPACK.")

if(WITH_IPOPT)
  find_package(IPOPT)
endif()
add_feature_info(ipopt-interface IPOPT_FOUND "Interface to the NLP solver Ipopt.")

if(WITH_KNITRO)
  find_package(KNITRO)
endif()
add_feature_info(knitro-interface KNITRO_FOUND "Interface to the NLP solver KNITRO.")

if(WITH_CPLEX)
  find_package(CPLEX)
endif()
add_feature_info(cplex-interface CPLEX_FOUND "Interface to the QP solver CPLEX.")

if(WITH_SNOPT)
  find_package(SNOPT)
endif()
add_feature_info(snopt-interface SNOPT_FOUND "Interface to SNOPT.")

# if(WITH_GSL)
#   find_package (GSL)
# endif()

# Use the included qpOASES
if(WITH_QPOASES AND BLAS_FOUND AND LAPACK_FOUND)
  set(QPOASES_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external_packages/qpOASES/include)
  set(QPOASES_LIBRARIES casadi_qpoases)
  set(QPOASES_FOUND true)
endif()
add_feature_info(qpoases-interface QPOASES_FOUND "Interface to the active-set QP solver qpOASES.")

if(WITH_OOQP AND BLAS_FOUND AND HSL_FOUND)
  find_package(OOQP)
endif()
add_feature_info(ooqp-interface OOQP_FOUND "Interface to the QP solver OOQP (requires BLAS and HSL libraries).")

if(Fortran_language_works AND WITH_SQIC)
  find_package(SQIC)
endif()
add_feature_info(sqic-interface SQIC_FOUND "Interface to the QP solver SQIC.")

if(Fortran_language_works AND WITH_SLICOT)
  find_package(SLICOT)
endif()
add_feature_info(slicot-interface SLICOT_FOUND "Interface to the controls library SLICOT.")


if(WITH_WORHP AND LIBXML2_FOUND  AND WSMP_FOUND AND METIS_FOUND AND BLAS_FOUND AND LAPACK_FOUND)
  find_package(WORHP)
endif()
add_feature_info(worhp-inteface WORHP_FOUND "Interface to the NLP solver Worhp (requires LibXml2, wsmp, metis, blas, lapack).")

# Use the included dsdp
if(WITH_DSDP AND BLAS_FOUND AND LAPACK_FOUND)
  set(DSDP_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external_packages/DSDP/include)
  set(DSDP_LIBRARIES casadi_dsdp)
  set(DSDP_FOUND true)
endif()
add_feature_info(dsdp-interface DSDP_FOUND "Interface to the interior point SDP solver DSDP (requires BLAS and LAPACK).")

if(WITH_PROFILING)
  find_library(RT rt)
  add_definitions(-DWITH_PROFILING)
endif()

add_subdirectory(external_packages)
add_subdirectory(casadi)
add_subdirectory(docs/examples)
add_subdirectory(experimental/greg EXCLUDE_FROM_ALL)
add_subdirectory(experimental/joel EXCLUDE_FROM_ALL)
#add_subdirectory(experimental/andrew EXCLUDE_FROM_ALL)
add_subdirectory(misc)

find_package(SWIG)
  if(SWIG_FOUND)
    add_subdirectory(swig EXCLUDE_FROM_ALL) # the swig interfaces are not compiled when typing "make"
    add_custom_target(python DEPENDS _python_modules)
    add_custom_target(json DEPENDS _json_modules)
endif()

if(WITH_DOC)
  set(DOXYFILE_LATEX "NO")
  set(DOXYFILE_IN_DIR ${CMAKE_SOURCE_DIR}/docs/api)
  set(DOXYFILE_OUTPUT_DIR  ${CMAKE_SOURCE_DIR}/docs/api)
  set(DOXYFILE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/symbolic)
  include(UseDoxygen OPTIONAL)
  add_subdirectory(docs/tutorials EXCLUDE_FROM_ALL)
  add_custom_target(tutor)
endif()

add_subdirectory(docs/api/examples/ctemplate)

# Print a summary
feature_summary(WHAT ALL)
feature_summary(WHAT ENABLED_FEATURES VAR feature_list)

string(REPLACE "\n" "\\n" feature_list "${feature_list}")

# Pass meta information
configure_file(
  ${CMAKE_SOURCE_DIR}/casadi/symbolic/casadi_meta.cpp.cmake
  ${CMAKE_SOURCE_DIR}/casadi/symbolic/casadi_meta.cpp
  ESCAPE_QUOTES
)

if(EXTRA_CMAKE)
  include(${EXTRA_CMAKE})
endif()
