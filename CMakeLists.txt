cmake_minimum_required(VERSION 2.8.12)
project(casadi C CXX)

set(CASADI_MAJOR_VERSION 3)
set(CASADI_MINOR_VERSION 5)
set(CASADI_PATCH_VERSION 5)
set(CASADI_IS_RELEASE 0)

set(CASADI_VERSION ${CASADI_MAJOR_VERSION}.${CASADI_MINOR_VERSION}.${CASADI_PATCH_VERSION})
if(NOT CASADI_IS_RELEASE)
  set(CASADI_VERSION "${CASADI_VERSION}+")
endif()
set(CASADI_VERSION "${CASADI_VERSION}")
set(PACKAGE_NAME "CasADi")
set(PACKAGE_VERSION ${CASADI_VERSION})
set(PACKAGE_SO_VERSION "3.6")
set(PACKAGE_DESCRIPTION "CasADi is a symbolic framework for automatic differentation and numeric optimization")
set(PACKAGE_AUTHOR "Joel Andersson and Joris Gillis")
set(PACKAGE_MAINTAINER "Joel Andersson and Joris Gillis")
set(PACKAGE_URL "http://casadi.org")

############################################################
####################### Policies ###########################
############################################################
# Enforce cmake2.4+ style behaviour for
# mixed filename/full path linker directives
cmake_policy(SET CMP0003 NEW)

# cmake_policy propagation rules
if(POLICY CMP0011)
cmake_policy(SET CMP0011 OLD)
endif()

# cmake_policy avoid escaping in add_definitions
cmake_policy(SET CMP0005 OLD)

# cmake_policy use MACOSX_RPATH=ON by default
if(POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW)
endif()

# cmake_policy how to interpret string arguments
if(POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif()

#############################################################
########################## Version ##########################
#############################################################
# For systems with git installed, find out revision and description.
execute_process(COMMAND git rev-parse HEAD
  OUTPUT_VARIABLE git_revision
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
execute_process(COMMAND git describe --first-parent HEAD
  OUTPUT_VARIABLE git_describe
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
string(REPLACE "+-" "+" git_describe "${git_describe}")
string(REPLACE "-g" "." git_describe "${git_describe}")
set(PACKAGE_VERSION_FULL "${git_describe}")

# overwrite version numbers for official release
if((PACKAGE_VERSION_FULL STREQUAL "") OR (NOT PACKAGE_VERSION MATCHES "\\+"))
  set(PACKAGE_VERSION_FULL "${PACKAGE_VERSION}")
endif()

# We use the FeatureSummary module to get a more readable CMake output
include(FeatureSummary)
set_property(GLOBAL APPEND PROPERTY FeatureSummary_PKG_TYPES BUILD)
if(NOT COMMAND add_feature_info)
  function(add_feature_info _name _enabled _description)
    message(STATUS "Detected: ${_name} ${_description}")
  endfunction()
  function(feature_summary)
  endfunction()
endif()

# Default build type is Release
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (default Release)" FORCE)
endif()

########################################################################
############################# Fortran ##################################
########################################################################
option(FORTRAN_REQUIRED "Make fortran a required dependency" OFF)

# http://public.kitware.com/Bug/print_bug_page.php?bug_id=9220
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(language_support)
workaround_9220(Fortran Fortran_language_works)
if(FORTRAN_REQUIRED)
  enable_language(Fortran)
  set(Fortran_language_works ON)
else()
  if(Fortran_language_works)
    enable_language(Fortran OPTIONAL)
  endif()
endif()

########################################################################
############################# C++11 ####################################
########################################################################
# Check if the compiler supports C++11
foreach(CXX11FLAG_TEST "" "-std=c++11" "-std=c++0x")
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
else()
  message(FATAL_ERROR "Cannot enable C++11 features.")
endif()

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

# Check if mkstemps is available
if(WIN32)
  add_definitions(-DHAVE_SIMPLE_MKSTEMPS)
else()
  file(WRITE ${CMAKE_BINARY_DIR}/mkstemps_test.cxx
      "#include <cstdlib>\n"
      "#include <unistd.h>\n"
      "int main() {\n"
      "  char tmp[] = \"mkstempstestXXXXXX.so\";\n"
      "  int flag = mkstemps(tmp, 3);\n"
      "  return flag;\n"
      "}\n")
  try_compile(HAVE_MKSTEMPS ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/mkstemps_test.cxx
              OUTPUT_VARIABLE HAVE_MKSTEMPS_OUTPUT)
  file(WRITE ${CMAKE_BINARY_DIR}/simple_mkstemps_test.cxx
      "#include <sys/stat.h>\n"
      "#include <fcntl.h>\n"
      "#include <errno.h>\n"
      "int main() {\n"
      "  int fd = open(\"test.c\", O_RDWR | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);\n"
      "}\n")
  try_compile(HAVE_SIMPLE_MKSTEMPS ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/simple_mkstemps_test.cxx
              OUTPUT_VARIABLE HAVE_SIMPLE_MKSTEMPS_OUTPUT)
  if(HAVE_MKSTEMPS)
    add_definitions(-DHAVE_MKSTEMPS)
  else()
    if(HAVE_SIMPLE_MKSTEMPS)
      add_definitions(-DHAVE_SIMPLE_MKSTEMPS)
    else()
      message(WARNING "mkstemp is not available, falling back on tmpnam. "
                      "Test compilation returned the following:\n"
                      "${HAVE_MKSTEMPS_OUTPUT}"
                      "Test compilation returned the following:\n"
                      "${HAVE_SIMPLE_MKSTEMPS_OUTPUT}")
   endif()
  endif()
endif()

############################################################
####################### macports python ####################
############################################################
# find the macports python libs if they exist
if(APPLE)
  set(CMAKE_LIBRARY_PATH /opt/local/lib ${CMAKE_LIBRARY_PATH})
endif()

############################################################
################# WITH_COMMON option #######################
############################################################
option(WITH_COMMON "Build common, stable interfaces, if found" OFF)
if(WITH_COMMON)
  # Look for LAPACK
  find_package(LAPACK QUIET)
  set(WITH_LAPACK_DEF ${LAPACK_FOUND})
  message("WITH_LAPACK_DEF: ${LAPACK_FOUND}")
  # qpOASES and blockSQP built if LAPACK available
  set(WITH_QPOASES_DEF ${LAPACK_FOUND})
  set(WITH_BLOCKSQP_DEF ${LAPACK_FOUND})
  set(WITH_SUPERSCS_DEF ${LAPACK_FOUND})
  # Look for IPOPT
  find_package(IPOPT QUIET)
  set(WITH_IPOPT_DEF ${IPOPT_FOUND})
else()
  # Disabled by default
  set(WITH_LAPACK_DEF OFF)
  set(WITH_QPOASES_DEF OFF)
  set(WITH_BLOCKSQP_DEF OFF)
  set(WITH_IPOPT_DEF OFF)
  set(WITH_SUPERSCS_DEF OFF)
endif()

set(DEPENDENT_LIBRARIES)


include(ExternalProject)


set(CROSS_FLAGS "")
set(CROSS_FLAGS_FORTRAN "")
set(CROSS_FLAGS_CMAKE)
if(CMAKE_CROSS_COMPILING)
  set(CROSS_FLAGS "--host=$ENV{CROSS_TRIPLE}")
  if("$ENV{CROSS_TRIPLE}" MATCHES "mingw32")
    # https://github.com/casadi/artifacts/issues/3
    set(CROSS_FLAGS_FORTRAN "lt_cv_deplibs_check_method=pass_all")
  endif()
  set(CROSS_FLAGS_CMAKE "-DCMAKE_TOOLCHAIN_FILE=$ENV{CMAKE_TOOLCHAIN_FILE}")
endif()

set(EXTERNAL_CMAKE_COMMAND "${CMAKE_COMMAND}")
if("$ENV{CROSS_TRIPLE}" MATCHES "mingw32")
  string(REPLACE "-gcc" "-cmake" EXTERNAL_CMAKE_COMMAND "$ENV{CC}")
endif()
message("EXTERNAL_CMAKE_COMMAND: ${EXTERNAL_CMAKE_COMMAND}")

set(WHOLE_ARCHIVE_ON "")
set(WHOLE_ARCHIVE_OFF "")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(WHOLE_ARCHIVE_ON "-Wl,--whole-archive ")
  set(WHOLE_ARCHIVE_OFF "-Wl,--no-whole-archive ")
  if(APPLE)
    set(WHOLE_ARCHIVE_ON "-Wl,-all_load ")
    set(WHOLE_ARCHIVE_OFF "-Wl,-noall_load")
  endif()
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(WHOLE_ARCHIVE_ON "-Wl,-force_load, ")
  set(WHOLE_ARCHIVE_OFF "")
  if(APPLE)
    set(WHOLE_ARCHIVE_ON "-Wl,-all_load ")
    set(WHOLE_ARCHIVE_OFF "-Wl,-noall_load ")
  endif()
endif()


message("CXX $ENV{CXX}")
message("FC $ENV{FC}")
message("CC $ENV{CC}")



#####################################################
######################### swig ######################
#####################################################
option(WITH_PYTHON "Compile the Python front-end" OFF)
option(WITH_MATLAB "Compile the MATLAB front-end" OFF)
option(WITH_OCTAVE "Compile the Octave front-end (experimental)" OFF)
option(WITH_JSON "Compile the JSON front-end" OFF)
option(WITH_OCTAVE_IMPORT "[DEPRECATED] Copy import.m file to install directory" OFF)
option(WITH_PYTHON3 "Use Python 3.*" OFF)

############################################################
################# Misc options #############################
############################################################
option(INSTALL_INTERNAL_HEADERS "Install all headers" OFF)
option(ENABLE_EXPORT_ALL "Export all symbols to a shared library" OFF)
option(ENABLE_STATIC "Compile static libraries (.a on Linux)" OFF)
option(ENABLE_SHARED "Compile shared libraries (.so on Linux)" ON)
option(SWIG_EXPORT "Export SWIG" OFF)
option(SWIG_IMPORT "Export SWIG" OFF)

# enabling OpenMP support if requested
option(WITH_OPENMP "Compile with parallelization support using OpenMP" OFF)
if(WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DWITH_OPENMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -DWITH_OPENMP")
  endif()
endif()

# enabling Threads support if requested
option(WITH_THREAD "Compile with parallelization support using POSIX Threads" OFF)
option(WITH_THREAD_MINGW "For MINGW, use mingw-std-threads" OFF)
if(WITH_THREAD)
  add_definitions(-DCASADI_WITH_THREAD)
endif()
if(MINGW AND WITH_THREAD_MINGW)
  add_definitions(-DCASADI_WITH_THREAD_MINGW)
else()
  if(WITH_THREAD)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  endif()
endif()


# OpenCL
option(WITH_OPENCL "Compile with OpenCL support (experimental)" OFF)
if(WITH_OPENCL)
  # Core depends on OpenCL for GPU calculations
  find_package(OpenCL REQUIRED)
  add_definitions(-DWITH_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIRS})
endif()
add_feature_info(opencl-support WITH_OPENCL "Enable just-in-time compiliation to CPUs and GPUs with OpenCL.")

# Enable: RTLD_DEEPBIND
option(WITH_DEEPBIND "Load plugins with RTLD_DEEPBIND (can be used to resolve conflicting libraries in e.g. MATLAB)" ON)
if(WITH_DEEPBIND)
  add_definitions(-DWITH_DEEPBIND)
endif()

# Suppress warnings for refcounting failures (at program closure)
option(WITH_REFCOUNT_WARNINGS "Issue warnings upon reference counting failure" OFF)
if(WITH_REFCOUNT_WARNINGS)
  add_definitions(-DWITH_REFCOUNT_WARNINGS)
endif()

# Have an so version?
option(WITH_SO_VERSION "Use an so version for the library (version suffix) when applicable" ON)

set(WITH_COPYSIGN_UNDEF_DEFAULT OFF)
if(WITH_PYTHON AND MSVC)
  set(WITH_COPYSIGN_UNDEF_DEFAULT OFF)
endif()
option(WITH_COPYSIGN_UNDEF "On Windows, in Python.h #defines copysign to be _copysign. This option inserts an #undef" WITH_COPYSIGN_UNDEF_DEFAULT)

############################################################
################# Installation location ####################
############################################################
# Flat hierarchy if self-contained and on Windows
option(WITH_SELFCONTAINED "Make the install directory self-contained" OFF)
if(WITH_SELFCONTAINED OR (WIN32 AND NOT CYGWIN))
  set(DEF_LIB_PREFIX casadi)
  set(DEF_BIN_PREFIX casadi)
  set(DEF_CMAKE_PREFIX casadi/cmake)
  set(DEF_PYTHON_PREFIX) # Note: Files installed in a casadi subdirectory
  set(DEF_MATLAB_PREFIX casadi)
  set(DEF_INCLUDE_PREFIX casadi/include)
else()
  set(DEF_LIB_PREFIX lib)
  set(DEF_BIN_PREFIX bin)
  set(DEF_CMAKE_PREFIX lib/cmake/casadi)
  set(DEF_PYTHON_PREFIX python)
  if (WITH_PYTHON)
    set(PYTHON_VERSION_REQ "")
    if (WITH_PYTHON3)
      set(PYTHON_VERSION_REQ "3")
    endif()
    find_package(PythonInterp ${PYTHON_VERSION_REQ} REQUIRED)
    execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c "from __future__ import print_function; from distutils.sysconfig import get_python_lib; print(get_python_lib())"
      OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(DEF_PYTHON_PREFIX ${PYTHON_SITE_PACKAGES})
  endif()
  set(DEF_MATLAB_PREFIX matlab)
  set(DEF_INCLUDE_PREFIX include)
endif()

# Installation directories are cache entries that can be modified
set(LIB_PREFIX ${DEF_LIB_PREFIX} CACHE PATH "Installation directory for libraries")
set(BIN_PREFIX ${DEF_BIN_PREFIX} CACHE PATH "Installation directory for executables")
set(INCLUDE_PREFIX ${DEF_INCLUDE_PREFIX} CACHE PATH "Installation directory for header files")
set(CMAKE_PREFIX ${DEF_CMAKE_PREFIX} CACHE PATH "Installation directory for CMake files")
set(PYTHON_PREFIX ${DEF_PYTHON_PREFIX} CACHE PATH "Installation directory for Python files")
set(MATLAB_PREFIX ${DEF_MATLAB_PREFIX} CACHE PATH "Installation directory for MATLAB/Octave files")


install(DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/bin/" OPTIONAL DESTINATION ${LIB_PREFIX})
install(DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/lib/" OPTIONAL DESTINATION ${LIB_PREFIX})
if(CMAKE_CROSS_COMPILING)
  install(DIRECTORY ${CMAKE_FIND_ROOT_PATH}/bin/ DESTINATION ${LIB_PREFIX} FILES_MATCHING PATTERN "*${CMAKE_SHARED_LIBRARY_SUFFIX}")
endif()

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE PYTHON MATLAB)
  set(var ${p}_PREFIX)
  if(${var})
    if(NOT IS_ABSOLUTE "${${var}}")
      # Turn relative path into absolute
      set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
    endif()
  else()
    # Flat installation hierarchy
    set(${var} "${CMAKE_INSTALL_PREFIX}")
  endif()
endforeach()

# Set rpath
if(WITH_SELFCONTAINED)
  if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path/")
  else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN")
  endif()
else()
  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # the RPATH to be used when installing, but only if it's not a system directory
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${LIB_PREFIX}" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${LIB_PREFIX}")
  endif()
endif()


################################################################################
############################### compiler warnings ##############################
################################################################################
# add compiler warnings
option(WITH_EXTRA_WARNINGS "Add extra warnings (like -Wall -Wextra)" OFF)
option(WITH_WERROR "Add -Werror" OFF)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if(WITH_EXTRA_WARNINGS)
    set(EXTRA_CASADI_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-sign-compare -Wno-attributes")
    set(EXTRA_CASADI_C_FLAGS "-Wall -Wextra -Wformat-security -Wno-unused-parameter -Wno-maybe-uninitialized -Wno-unused-but-set-variable -Wno-unused-function -Wno-unused-value")
  endif()
  if(WITH_WERROR)
    set(EXTRA_CASADI_CXX_FLAGS "${EXTRA_CASADI_CXX_FLAGS} -Werror")
    set(EXTRA_CASADI_C_FLAGS "${EXTRA_CASADI_C_FLAGS} -Werror")
  endif()

elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  if(WITH_EXTRA_WARNINGS)
    set(EXTRA_CASADI_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-sign-compare -Wno-mismatched-tags -Wno-tautological-compare")
    set(EXTRA_CASADI_C_FLAGS "-Wall -Wextra -Wformat-security -Wno-unused-parameter -Wno-unused-function -Wno-unused-value")
  endif()
  if(WITH_WERROR)
    set(EXTRA_CASADI_CXX_FLAGS "${EXTRA_CASADI_CXX_FLAGS} -Werror")
    set(EXTRA_CASADI_C_FLAGS "${EXTRA_CASADI_C_FLAGS} -Werror")
  endif()

elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  if(WITH_EXTRA_WARNINGS)
    set(EXTRA_CASADI_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-sign-compare")
    set(EXTRA_CASADI_C_FLAGS "-Wall -Wextra -Wformat-security -Wno-unused-parameter -Wno-maybe-uninitialized -Wno-unused-but-set-variable -Wno-unused-function -Wno-unused-value")
  endif()
  if(WITH_WERROR)
    set(EXTRA_CASADI_CXX_FLAGS "${EXTRA_CASADI_CXX_FLAGS} -Werror")
    set(EXTRA_CASADI_C_FLAGS "${EXTRA_CASADI_C_FLAGS} -Werror")
  endif()

elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # I don't know the MSVC warnings, so do nothing for now
endif()
#######################################################################
############################### Extra runtime checks ##################
#######################################################################
option(WITH_EXTRA_CHECKS "Compile in extra runtime checks (for developers)" OFF)
if(WITH_EXTRA_CHECKS)
   add_definitions(-DWITH_EXTRA_CHECKS)
endif()


#######################################################################
############################# -fPIC ###################################
#######################################################################
# the following commands are needed to fix a problem with the libraries
# for linux 64 bits
if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64" OR "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64")
  message(STATUS "x86_64 architecture detected - setting flag -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
endif()

#######################################################################
########################### code coverage #############################
#######################################################################
option(WITH_COVERAGE "Create coverage report" OFF)
if(WITH_COVERAGE)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fprofile-arcs -ftest-coverage")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage")
set(CMAKE_EXE_LINKER_FLAGS  " -lgcov -fprofile-arcs --coverage ${CMAKE_EXE_LINKER_FLAGS}")
endif()

if(MINGW)
  # Circumventing a bug in MinGW g++ v4.7.2, evoked by 752fa89355ffa
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-ipa-cp-clone")
endif()

# Tentative fix for https://github.com/casadi/casadi/issues/1665
if(MINGW OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel"))
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fkeep-inline-functions")
endif()

# Plugin version
add_definitions(-DCASADI_VERSION=31)

# Necessary for Visual C++
add_definitions(-D_USE_MATH_DEFINES)
add_definitions(-D_SCL_SECURE_NO_WARNINGS)

# Dynamic loading (needed for loading plugins and external functions)
option(WITH_DL "Enable dynamic loading of functions" ON)
if(WITH_DL)
  add_definitions(-DWITH_DL)
endif()
add_feature_info(dynamic-loading WITH_DL "Compile with support for dynamic loading of generated functions (needed for ExternalFunction)")

# Dynamically loading FMU libraries
option(WITH_FMU "Enable dynamic loading of functions" ON)
if(WITH_FMU)
  add_definitions(-DWITH_FMU)
  set(FMU_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external_packages/FMI-Standard-2.0.2/headers)
endif()
add_feature_info(dynamic-loading WITH_FMU "Compile with support for dynamic loading of FMU libraries")

# Include support for deprecated features (to be removed in the next release)
option(WITH_DEPRECATED_FEATURES "Compile with syntax that is scheduled to be deprecated" ON)
if (WITH_DEPRECATED_FEATURES)
  add_definitions(-DWITH_DEPRECATED_FEATURES)
endif()

include_directories(.)
include_directories(${PROJECT_BINARY_DIR})

#######################################################################
################# 3rd party interface #################################
#######################################################################

option(WITH_BUILD_REQUIRED "Build any requirements that are not found on your system")

# Sundials: Solve initial-value problems in ODE and DAE
option(WITH_BUILD_SUNDIALS "Compile the included source code for Sundials 2.6.1" ON)
option(WITH_SUNDIALS "Compile the interface to Sundials" ON)
if(WITH_SUNDIALS)
  if(WITH_BUILD_SUNDIALS)
    # build the included sundials
    set(SUNDIALS_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external_packages/casadi-sundials)
    set(SUNDIALS_LIBRARIES casadi_sundials)
  else()
    # try to find system sundials
    add_definitions(-DWITH_SYSTEM_SUNDIALS)
    find_package(SUNDIALS REQUIRED)
  endif()
endif()
add_feature_info(sundials-interface WITH_SUNDIALS "Interface to the ODE/DAE integrator suite SUNDIALS.")

# CSparse: A sparse direct linear solver (included in the source)
option(WITH_BUILD_CSPARSE "Compile the included source code for CSparse" ON)
option(WITH_CSPARSE "Compile the interface to CSparse" ON)
if(WITH_CSPARSE)
  if(WITH_BUILD_CSPARSE)
    # build the included csparse
    set(CSPARSE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external_packages/CSparse/Include)
    set(CSPARSE_LIBRARIES casadi_csparse)
  else()
    # try to find system csparse
    find_package(CSPARSE REQUIRED)
  endif()
endif()
add_feature_info(csparse-interface WITH_CSPARSE "Interface to the sparse direct linear solver CSparse.")

# BLASFEO
option(WITH_BLASFEO "Compile the interface to BLASFEO" OFF)
set(BUILD_BLASFEO_VERSION "edf92b396adddd9e548b9786f87ad290a0971329")
option(WITH_BUILD_BLASFEO "Build BLASFEO (BUILD_BLASFEO_VERSION=${BUILD_BLASFEO_VERSION}) from downloaded source." OFF)

# HPIPM
option(WITH_HPIPM "Compile the interface to HPIPM" OFF)
set(BUILD_HPIPM_VERSION "b140f3633ef383cdf481ba4e8c94bdc21596c47b")
option(WITH_BUILD_HPIPM "Build HPIPM (BUILD_HPIPM_VERSION=${BUILD_HPIPM_VERSION}) from downloaded source." OFF)
set(BUILD_HPIPM_DEPENDENCIES BLASFEO)

# SuperSCS: An conic solver
option(WITH_SUPERSCS "Compile the interface to SuperSCS (the source code for SuperSCS is included)" ${WITH_SUPERSCS_DEF})
add_feature_info(superscs-interface WITH_SUPERSCS "Interface to Conic solver SuperSCS.")

# OSQP: An operator-splitting QP solver
option(WITH_OSQP "Compile the interface to OSQP (the source code for OSQP is included)" OFF)
add_feature_info(osqp-interface WITH_OSQP "Interface to QP solver OSQP.")

# PROXQP: A QP solver based on proximal methods
option(WITH_PROXQP "Compile the interface to PROXQP (the source code for PROXQP is included)" OFF)
add_feature_info(proxqp-interface WITH_PROXQP "Interface to QP solver PROXQP.")

# TinyXML: An XML reader, needed e.g. by the Modelica interface
option(WITH_BUILD_TINYXML "Compile the included TinyXML source code" ON)
option(WITH_TINYXML "Compile the interface to TinyXML" ON)
if(WITH_TINYXML)
  if(WITH_BUILD_TINYXML)
    # build the included tinyxml
    set(TINYXML_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external_packages/tinyxml-2.6.2)
    set(TINYXML_LIBRARIES casadi_tinyxml)
  else()
    # try to find system tinyxml
    find_package(TINYXML REQUIRED)
  endif()
endif()
add_feature_info(tinyxml-interface WITH_TINYXML "Interface to the XML parser TinyXML.")

# DSDP
option(WITH_BUILD_DSDP "Compile the the included source code for DSDP" OFF)
option(WITH_DSDP "Compile the interface to DSDP" OFF)
if(WITH_DSDP)
  if(WITH_BUILD_DSDP)
    # build the included dsdp
    set(DSDP_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external_packages/DSDP/include)
    set(DSDP_LIBRARIES casadi_dsdp)
  else()
    # try to find system dsdp
    find_package(DSDP REQUIRED)
  endif()
endif()
add_feature_info(dsdp-interface WITH_DSDP "Interface to the interior point SDP solver DSDP.")

# CLANG: A just-in-time compiler
option(OLD_LLVM "Use pre-3.5 LLVM" OFF)
option(WITH_CLANG "Compile the interface to clang JIT" OFF)
if(WITH_CLANG)
  # try to find system CLANG
  find_package(CLANG REQUIRED)
endif()
add_feature_info(clang-interface WITH_CLANG "Interface to the Clang JIT compiler.")

# Lapack: Dense linear solvers
option(WITH_LAPACK "Compile the interface to LAPACK" ${WITH_LAPACK_DEF})
option(WITH_BUILD_LAPACK "Download and install OpenBLAS for LAPACK+BLAS" OFF)

# qpOASES: An active-set QP solver
option(WITH_QPOASES "Compile the interface to qpOASES (the source code for qpOASES is included)" ${WITH_QPOASES_DEF})
option(WITH_NO_QPOASES_BANNER "Add -D__NO_COPYRIGHT__ to qpOASES definitions" OFF)
if(WITH_QPOASES)
  if (NOT WITH_LAPACK)
     message(FATAL_ERROR "WITH_QPOASES requires WITH_LAPACK.")
  endif()
  set(QPOASES_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external_packages/qpOASES/include)
  set(QPOASES_LIBRARIES casadi_qpoases)
endif()
add_feature_info(qpoases-interface WITH_QPOASES "Interface to the active-set QP solver qpOASES.")

# blockSQP: A block-structure-exploiting SQP method
option(WITH_BLOCKSQP "Compile the interface to blockSQP (the source code for blockSQP is included)" ${WITH_BLOCKSQP_DEF})
if(WITH_BLOCKSQP)
  if (NOT WITH_QPOASES)
    message(FATAL_ERROR "WITH_BLOCKSQP requires WITH_QPOASES.")
  endif()
endif()
add_feature_info(blocksqp-interface WITH_BLOCKSQP "Interface to the NLP solver blockSQP.")

# IPOPT: A sparse interior point NLP solver with filter line-search
option(WITH_IPOPT "Compile the interface to IPOPT" ${WITH_IPOPT_DEF})
set(BUILD_IPOPT_VERSION "3.14.4")
option(WITH_BUILD_IPOPT "Build IPOPT (BUILD_IPOPT_VERSION=${BUILD_IPOPT_VERSION}) from downloaded source." OFF)
set(BUILD_IPOPT_DEPENDENCIES MUMPS)

# Knitro: An NLP solver
option(WITH_KNITRO "Compile the interface to KNITRO" OFF)
if(WITH_KNITRO)
  find_package(KNITRO REQUIRED)
endif()
add_feature_info(knitro-interface WITH_KNITRO "Interface to the NLP solver KNITRO.")

# CPLEX: A sparse mixed-integer QP solver
option(WITH_CPLEX "Compile the interface to CPLEX" OFF)
option(WITH_CPLEX_SHARED "Use shared cplex libraries" OFF)
if(WITH_CPLEX)
  find_package(CPLEX REQUIRED)
endif()
add_feature_info(cplex-interface WITH_CPLEX "Interface to the QP solver CPLEX.")

# Bonmin: An MINLP solver
option(WITH_BONMIN "Compile the interface to BONMIN" OFF)
set(BUILD_BONMIN_VERSION "1.8.8")
option(WITH_BUILD_BONMIN "Build BONMIN (BUILD_BONMIN_VERSION=${BUILD_BONMIN_VERSION}) from downloaded source." OFF)
set(BUILD_BONMIN_DEPENDENCIES IPOPT CBC)

# CBC: A sparse mixed-integer linear solver
option(WITH_CBC "Compile the CBC interface" OFF)
set(BUILD_CBC_VERSION "2.10.6")
option(WITH_BUILD_CBC "Build CBC (BUILD_CBC_VERSION=${BUILD_CBC_VERSION}) from downloaded source." OFF)
set(BUILD_CBC_DEPENDENCIES CLP)

# CLP: A sparse linear solver
option(WITH_CLP "Compile the CLP interface" OFF)
set(BUILD_CLP_VERSION "1.17.7")
set(BUILD_COINUTILS_VERSION "2.11.6")
set(BUILD_CGL_VERSION "0.60.4")
set(BUILD_OSI_VERSION "0.108.7")
option(WITH_BUILD_CLP "Build CLP (BUILD_CLP_VERSION=${BUILD_CLP_VERSION}) from downloaded source." OFF)
set(BUILD_CLP_DEPENDENCIES MUMPS LAPACK)

# MUMPS: Sparse direct linear solver
option(WITH_MUMPS "Enable MUMPS interface" OFF)
set(BUILD_MUMPS_TP_VERSION "3.0.2")
set(BUILD_METIS_TP_VERSION "2.0.0")
option(WITH_BUILD_MUMPS "Build MUMPS (BUILD_MUMPS_TP_VERSION=${BUILD_MUMPS_TP_VERSION}) from downloaded source." OFF)
set(BUILD_MUMPS_DEPENDENCIES LAPACK)

# HSL: Sparse direct linear solvers
option(WITH_HSL "Enable HSL interface" OFF)
if(WITH_HSL)
  find_package(HSL REQUIRED)
endif()
add_feature_info(hsl-interface WITH_HSL "Interface to HSL.")


function(libraries2flags libraries out)
# libtool cannot handle absoulte paths (i.e. it omits the library from the list of dependencies)
  set(local "")
  foreach(LIB ${libraries})
    if (LIB MATCHES "^-")
      # Pass -l, -L as-is
      set(local "${local} ${LIB}")
    else()
      get_filename_component(LIB_DIR "${LIB}" DIRECTORY)
      get_filename_component(LIB_NAME "${LIB}" NAME)
      string(REGEX REPLACE "${CMAKE_SHARED_LIBRARY_PREFIX}([A-Za-z0-9_]+)${CMAKE_SHARED_LIBRARY_SUFFIX}" "\\1" LIB_NAME_MOD ${LIB_NAME})
      if(LIB_NAME STREQUAL LIB_NAME_MOD)
        set(local "${local} ${LIB}")
      else()
        set(local "${local} -L${LIB_DIR} -l${LIB_NAME_MOD}")
      endif()
    endif()
    set(${out} "${local}" PARENT_SCOPE)
  endforeach()
endfunction()

# Packages must be sorted from high to low in the dependency graph
foreach(PKG BONMIN IPOPT CBC CLP MUMPS HPIPM BLASFEO LAPACK)
  string(TOLOWER ${PKG} PKG_LOWER)
  if(WITH_${PKG})
    if(NOT WITH_BUILD_${PKG})
      get_property(PACKAGES_NOT_FOUND_BACKUP GLOBAL PROPERTY PACKAGES_NOT_FOUND)
      get_property(PACKAGES_FOUND_BACKUP GLOBAL PROPERTY PACKAGES_FOUND)
      find_package(${PKG})
    endif()
    if (WITH_BUILD_${PKG} OR (WITH_BUILD_REQUIRED AND NOT ${PKG}_FOUND))
      set(WITH_BUILD_${PKG} ON)

      # Add it to FeatureSummary
      get_property(PKG_HELP CACHE WITH_BUILD_${PKG} PROPERTY HELPSTRING)
      add_feature_info(${PKG_LOWER}-sourcebuild ON ${PKG_HELP})

      # Remove from PACKAGES_NOT_FOUND
      set_property(GLOBAL PROPERTY PACKAGES_NOT_FOUND ${PACKAGES_NOT_FOUND_BACKUP})
      set_property(GLOBAL PROPERTY PACKAGES_FOUND ${PACKAGES_FOUND_BACKUP})

      # We are going to need build dependencies
      foreach(DEP ${BUILD_${PKG}_DEPENDENCIES})
        set(WITH_${DEP} ON)
      endforeach()
      set(${PKG}_FOUND ON)
    endif()
    if(NOT ${PKG}_FOUND)
      message(FATAL_ERROR "${PKG} not found on your system yet required. Set WITH_BUILD_${PKG} or WITH_BUILD_REQUIRED ON to build this required package from source.")
    endif()
  endif()
endforeach()

if(WITH_LAPACK)
  if (WITH_BUILD_LAPACK)
    ExternalProject_Add(openblas-external
      GIT_REPOSITORY https://github.com/jgillis/OpenBLAS.git
      GIT_TAG ee15f53c9ef7148b05e1b2ff80003df1685f0f10
      PREFIX "${CMAKE_BINARY_DIR}/external_projects"
      UPDATE_COMMAND ""
      CMAKE_COMMAND "${EXTERNAL_CMAKE_COMMAND}"
      CMAKE_ARGS -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_LIBDIR=lib -DTARGET=CORE2 -DUSE_THREAD=ON -DNUM_THREADS=16 -DBUILD_SINGLE=ON -DBUILD_DOUBLE=ON -DBUILD_COMPLEX=ON -DBUILD_COMPLEX16=ON -DDYNAMIC_ARCH=OFF -DUSE_OPENMP=OFF -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
    set(LAPACK_LIBRARIES "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}openblas${CMAKE_SHARED_LIBRARY_SUFFIX}")
    set(LAPACK_DEPENDS_TARGET openblas-external)

    add_library(lapack SHARED IMPORTED)
    set_target_properties(lapack PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}openblas${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}openblas${CMAKE_IMPORT_LIBRARY_SUFFIX}"
    )
  else()
    add_library(lapack INTERFACE)
    target_link_libraries(lapack INTERFACE ${LAPACK_LIBRARIES})
  endif()
endif()
if(LAPACK_FOUND)
  libraries2flags("${LAPACK_LIBRARIES}" LAPACK_LFLAGS)
  message("LAPACK_LFLAGS: ${LAPBLASFEO_PATHACK_LFLAGS}")
endif()


if(WITH_BLASFEO)
  if (WITH_BUILD_BLASFEO)
    ExternalProject_Add(blasfeo-external
      GIT_REPOSITORY https://github.com/giaf/blasfeo.git
      GIT_TAG ${BUILD_BLASFEO_VBLASFEO_PATHERSION}
      PREFIX "${CMAKE_BINARY_DIR}/external_projects"
      UPDATE_COMMAND ""
      CMAKE_COMMAND "${EXTERNAL_CMAKE_COMMAND}"
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
    set(BLASFEO_PATH "${CMAKE_BINARY_DIR}/external_projects")
    set(BLASFEO_DEPENDS_TARGET blasfeo-external)

    add_library(blasfeo SHARED IMPORTED)
    set_target_properties(blasfeo PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/include"
    )
  else()
    #add_library(lapack INTERFACE)
    #target_link_libraries(lapack INTERFACE ${LAPACK_LIBRARIES})
  endif()
endif()

if(WITH_HPIPM)
  if (WITH_BUILD_HPIPM)
    ExternalProject_Add(hpipm-external
      DEPENDS ${BLASFEO_DEPENDS_TARGET}
      GIT_REPOSITORY https://github.com/jgillis/hpipm.git
      GIT_TAG ${BUILD_HPIPM_VERSION}
      PREFIX "${CMAKE_BINARY_DIR}/external_projects"
      UPDATE_COMMAND ""
      CMAKE_COMMAND "${EXTERNAL_CMAKE_COMMAND}"
      CMAKE_ARGS -DBLASFEO_PATH=${BLASFEO_PATH} -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
    set(LAPACK_LIBRARIES "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}openblas${CMAKE_SHARED_LIBRARY_SUFFIX}")
    set(LAPACK_DEPENDS_TARGET openblas-external)

    add_library(hpipm SHARED IMPORTED)
    target_link_libraries(hpipm INTERFACE blasfeo)
    set_target_properties(hpipm PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}openblas${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}openblas${CMAKE_IMPORT_LIBRARY_SUFFIX}"
    )
  else()
    #add_library(lapack INTERFACE)
    #target_link_libraries(lapack INTERFACE ${LAPACK_LIBRARIES})
  endif()
endif()


add_feature_info(lapack-interface WITH_LAPACK "Interface to LAPACK.")

if(WITH_MUMPS)
  if(WITH_BUILD_MUMPS)
    add_feature_info(metis-sourcebuild ON "Build METIS (BUILD_METIS_TP_VERSION=${BUILD_METIS_TP_VERSION}) from downloaded source.")
    ExternalProject_Add(metis-external
        GIT_REPOSITORY https://github.com/coin-or-tools/ThirdParty-Metis.git
        GIT_TAG releases/${BUILD_METIS_TP_VERSION}
        DEPENDS ${LAPACK_DEPENDS_TARGET}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS} --prefix=<INSTALL_DIR>

          # GKFree is not defined (variadic definition in utils.c); fails on MacOS target
          "CFLAGS=-Wno-implicit-function-declaration"
        PATCH_COMMAND COMMAND cd <SOURCE_DIR> && ./get.Metis)
    set(EXTERNAL_METIS_LIBS coinmetis)

    set(DEPENDENT_LIBRARIES ${DEPENDENT_LIBRARIES} ${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}${EXTERNAL_METIS_LIBS}${CMAKE_SHARED_LIBRARY_SUFFIX})

    ExternalProject_Add(mumps-external
        DEPENDS metis-external
        GIT_REPOSITORY https://github.com/coin-or-tools/ThirdParty-Mumps.git
        GIT_TAG releases/${BUILD_MUMPS_TP_VERSION}
        DEPENDS ${LAPACK_DEPENDS_TARGET} metis-external
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS}  --prefix=<INSTALL_DIR> "--with-metis-cflags=-I${CMAKE_BINARY_DIR}/external_projects/include/coin-or/metis" "--with-metis-lflags=-L${CMAKE_BINARY_DIR}/external_projects/lib -l${EXTERNAL_METIS_LIBS} -lm" "--with-lapack-lflags=${LAPACK_LFLAGS}" ${CROSS_FLAGS_FORTRAN}
        PATCH_COMMAND COMMAND cd <SOURCE_DIR> && ./get.Mumps)

    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/include/coin-or")
    add_library(mumps SHARED IMPORTED)
    set_target_properties(mumps PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}coinmumps${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}coinmumps${CMAKE_IMPORT_LIBRARY_SUFFIX}"
      INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/include/coin-or"
    )
  endif()
endif()
add_feature_info(mumps-interface WITH_MUMPS "Interface to MUMPS.")

if(WITH_CLP)
  if(WITH_BUILD_CLP)
    add_feature_info(coinutils-sourcebuild ON "Build COINUTILS (BUILD_COINUTILS_VERSION=${BUILD_COINUTILS_VERSION}) from downloaded source.")
    ExternalProject_Add(coinutils-external
          DEPENDS ${LAPACK_DEPENDS_TARGET}
          GIT_REPOSITORY https://github.com/coin-or/CoinUtils.git
          GIT_TAG releases/${BUILD_COINUTILS_VERSION}
          PREFIX "${CMAKE_BINARY_DIR}/external_projects"
          UPDATE_COMMAND ""
          CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS} "--with-lapack-lflags=${LAPACK_LFLAGS}" --prefix=<INSTALL_DIR> )
    add_feature_info(osi-sourcebuild ON "Build OSI (BUILD_OSI_VERSION=${BUILD_OSI_VERSION}) from downloaded source.")
    ExternalProject_Add(osi-external
        DEPENDS coinutils-external
        GIT_REPOSITORY https://github.com/coin-or/Osi.git
        GIT_TAG releases/${BUILD_OSI_VERSION}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS} --prefix=<INSTALL_DIR> ${CROSS_FLAGS_FORTRAN}
	# For unknown reasons, Coin-or/configure disables pkg-config by default for clang compilers 
	--enable-pkg-config
    )
    add_library(osi SHARED IMPORTED)
    set_target_properties(osi PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}Osi${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}Osi${CMAKE_IMPORT_LIBRARY_SUFFIX}"
    )
    ExternalProject_Add(clp-external
        DEPENDS osi-external mumps-external coinutils-external ${LAPACK_DEPENDS_TARGET}
        GIT_REPOSITORY https://github.com/coin-or/Clp.git
        GIT_TAG releases/${BUILD_CLP_VERSION}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS} "--with-lapack-lflags=${LAPACK_LFLAGS}" --prefix=<INSTALL_DIR>  ${CROSS_FLAGS_FORTRAN}
	# For unknown reasons, Coin-or/configure disables pkg-config by default for clang compilers 
	--enable-pkg-config
	)

    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/include/coin")
    add_library(clp-raw SHARED IMPORTED)
    set_target_properties(clp-raw PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIXBLASFEO_PATH}Clp${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}Clp${CMAKE_IMPORT_LIBRARY_SUFFIX}"
      INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/include/coin"
    )
    add_library(osi-clp SHARED IMPORTED)
    set_target_properties(osi-clp PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}OsiClp${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}OsiClp${CMAKE_IMPORT_LIBRARY_SUFFIX}"
    )
    add_library(clp INTERFACE)
    target_link_libraries(clp INTERFACE clp-raw osi-clp coinutils)
  endif()
endif()
add_feature_info(clp-interface WITH_CLP "Interface to the LP solver CLP.")


if(WITH_CBC)
  if(WITH_BUILD_CBC)
    add_feature_info(cgl-sourcebuild ON "Build CGL (BUILD_CGL_VERSION=${BUILD_CGL_VERSION}) from downloaded source.")
    ExternalProject_Add(cgl-external
        DEPENDS osi-external clp-external coinutils-external
        GIT_REPOSITORY https://github.com/coin-or/Cgl.git
        GIT_TAG releases/${BUILD_CGL_VERSION}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS} --prefix=<INSTALL_DIR> ${CROSS_FLAGS_FORTRAN}
	# For unknown reasons, Coin-or/configure disables pkg-config by default for clang compilers 
	--enable-pkg-config
	)
    add_library(coinutils SHARED IMPORTED)
    set_target_properties(coinutils PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}CoinUtils${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}CoinUtils${CMAKE_IMPORT_LIBRARY_SUFFIX}"
    )
    ExternalProject_Add(cbc-external
        DEPENDS coinutils-external cgl-external osi-external mumps-external clp-external ${LAPACK_DEPENDS_TARGET}
        GIT_REPOSITORY https://github.com/coin-or/Cbc.git
        GIT_TAG releases/${BUILD_CBC_VERSION}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS} --prefix=<INSTALL_DIR> "--with-lapack-lflags=${LAPACK_LFLAGS}"  ${CROSS_FLAGS_FORTRAN}
	# For unknown reasons, Coin-or/configure disables pkg-config by default for clang compilers 
	--enable-pkg-config
	)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/include/coin")
    add_library(cbc-raw SHARED IMPORTED)
    set_target_properties(cbc-raw PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}Cbc${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}Cbc${CMAKE_IMPORT_LIBRARY_SUFFIX}"
      INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/include/coin"
    )
    add_library(cbc-solver SHARED IMPORTED)
    set_target_properties(cbc-solver PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}CbcSolver${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}CbcSolver${CMAKE_IMPORT_LIBRARY_SUFFIX}"
    )
    add_library(cbc INTERFACE)
    target_link_libraries(cbc INTERFACE cbc-raw cbc-solver coinutils osi osi-clp)
  endif()
endif()
add_feature_info(cbc-interface WITH_CBC "Interface to the LP solver CBC.")

if(WITH_IPOPT)
  if(WITH_BUILD_IPOPT)
    ExternalProject_Add(ipopt-external
        DEPENDS ${LAPACK_DEPENDS_TARGET} mumps-external
        URL https://www.coin-or.org/download/source/Ipopt/Ipopt-${BUILD_IPOPT_VERSION}.zip
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS} --prefix=<INSTALL_DIR> "--with-mumps-cflags=-I${CMAKE_BINARY_DIR}/external_projects/include/coin-or/mumps" "--with-mumps-lflags=-L${CMAKE_BINARY_DIR}/external_projects/lib ${WHOLE_ARCHIVE_ON}-lcoinmumps ${WHOLE_ARCHIVE_OFF}" "--with-lapack-lflags= ${LAPACK_LFLAGS}" ${CROSS_FLAGS_FORTRAN}
	# For unknown reasons, Coin-or/configure disables pkg-config by default for clang compilers 
	--enable-pkg-config
	)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/include/coin-or")

    add_library(ipopt SHARED IMPORTED)
    set_target_properties(ipopt PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}ipopt${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}ipopt${CMAKE_IMPORT_LIBRARY_SUFFIX}"
      INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/include/coin-or"
    )
  endif()
endif()
add_feature_info(ipopt-interface WITH_IPOPT "Interface to the NLP solver Ipopt.")

if(WITH_BONMIN)
  if(WITH_BUILD_BONMIN)
    ExternalProject_Add(bonmin-external
        DEPENDS mumps-external metis-external coinutils-external osi-external clp-external cgl-external cbc-external ipopt-external
        GIT_REPOSITORY https://github.com/coin-or/Bonmin.git
        GIT_TAG releases/${BUILD_BONMIN_VERSION}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS} --prefix=<INSTALL_DIR>  ${CROSS_FLAGS_FORTRAN}
	# For unknown reasons, Coin-or/configure disables pkg-config by default for clang compilers 
	--enable-pkg-config
	)

    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/include/coin")
    add_library(bonmin-raw SHARED IMPORTED)
    set_target_properties(bonmin-raw PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}bonmin${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}bonmin${CMAKE_IMPORT_LIBRARY_SUFFIX}"
      INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/include/coin"
    )
    add_library(bonmin INTERFACE)
    target_link_libraries(bonmin INTERFACE bonmin-raw ipopt coinutils)
  endif()
endif()
add_feature_info(bonmin-interface WITH_BONMIN "Interface to the MINLP framework Bonmin.")


# SNOPT: An SQP method
option(WITH_SNOPT "Compile the SNOPT interface if SNOPT libraries can be found" OFF)
option(WITH_SNOPT_FORCE_DUMMY "Compile the SNOPT interface against a dummy library (for Debian distribution; breaks normal SNOPT build)" OFF)
option(WITH_BUILD_SNOPT_INTERFACE "Use the SNOPT C-interface shipped with casadi" ON)
if(WITH_SNOPT)
  if(NOT WITH_SNOPT_FORCE_DUMMY)
    find_package(SNOPT REQUIRED)
  endif()
  if(WITH_BUILD_SNOPT_INTERFACE)
    set(SNOPT_INTERFACE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external_packages/snopt-interface/include)
    set(SNOPT_INTERFACE_LIBRARIES casadi_snopt_interface)
  else()
    find_package(SNOPT_INTERFACE REQUIRED)
  endif()
endif()
add_feature_info(snopt-interface WITH_SNOPT "Interface to the NLP solver SNOPT.")

# HIGHS: A MILP / QP solver
option(WITH_HIGHS "Compile the HIGHS interface" OFF)
if(WITH_HIGHS)
  find_package(HIGHS REQUIRED)
endif()
add_feature_info(highs-interface WITH_HIGHS "Interface to the MILP / QP solver HIGHS.")

# OOQP: A sparse interior point QP solver
option(WITH_OOQP "Enable OOQP interface" OFF)
if(WITH_OOQP)
  find_package(OOQP REQUIRED)
endif()
add_feature_info(ooqp-interface WITH_OOQP "Interface to the QP solver OOQP (requires BLAS and HSL libraries).")

# GUROBI: A sparse mixed-integer QP solver
option(WITH_GUROBI "Compile the interface to GUROBI" OFF)
if(WITH_GUROBI)
  find_package(GUROBI REQUIRED)
endif()
add_feature_info(gurobi-interface WITH_GUROBI "Interface to the (mixed-integer) QP solver GUROBI")

# SQIC: A QP solver
option(WITH_SQIC "Enable SQIC interface" OFF)
if(WITH_SQIC)
  if (NOT Fortran_language_works)
    message(FATAL_ERROR "SQIC requires Fortran support.")
  endif()
  find_package(SQIC REQUIRED)
endif()
add_feature_info(sqic-interface WITH_SQIC "Interface to the QP solver SQIC.")

# AMPL: Interface to the AMPL solver library and interfaced solvers
option(WITH_AMPL "Compile the interface to AMPL" OFF)
add_feature_info(ampl-interface WITH_AMPL "Interface to the AMPL solver library.")

# Slicot: used for expm, dple
option(WITH_SLICOT "Enable SLICOT interface" OFF)
if(WITH_SLICOT)
  if (NOT Fortran_language_works)
    message(FATAL_ERROR "SLICOT requires Fortran support.")
  endif()
  if (NOT WITH_LAPACK)
     message(FATAL_ERROR "WITH_SLICOT requires WITH_LAPACK.")
  endif()
  find_package(SLICOT REQUIRED)
endif()
add_feature_info(slicot-interface WITH_SLICOT "Interface to the controls library SLICOT.")

# Worhp: A sparse SQP method
option(WITH_WORHP "Compile the WORHP interface" OFF)
if(WITH_WORHP)
  find_package(WORHP REQUIRED)
endif()
add_feature_info(worhp-interface WITH_WORHP "Interface to the NLP solver Worhp (requires gfortran, gomp).")

######################################################
##################### paths ##########################
######################################################
if(WIN32)
  # On Windows, build everything in the same directory
  set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
  set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})
else()
  # On non-Windows, use separate directories for build
  set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
  set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
endif()

#####################################################
###################### swig #########################
#####################################################
option(WITH_EXTENDING_CASADI "Compile a demonstration that shows how a project that depends on CasADi can be implemented." OFF)
if(WITH_PYTHON OR WITH_MATLAB OR WITH_OCTAVE OR WITH_JSON)
  add_subdirectory(swig)
endif()

option(WITH_LINT "Compile with support for linting (syntax checking)" OFF)
if(WITH_LINT)
  include(CppLint REQUIRED)
endif()

option(WITH_SPELL "Compile with support for spell-checking" OFF)
if(WITH_SPELL)
  include(SpellCheck REQUIRED)
endif()

option(WITH_CLANG_TIDY "Compile with support for clang-tidy" OFF)
if(WITH_CLANG_TIDY)
  include(ClangTidy REQUIRED)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  if(WITH_EXTRA_WARNINGS)
    # -Wall leaks into clang-tidy (through compile_commands.json)
    # silently overruling the 'Checks' field.
    # These warnings come from external_packages
    set(EXTRA_CASADI_CXX_FLAGS "${EXTRA_CASADI_CXX_FLAGS} -Wno-overloaded-virtual -Wno-deprecated-register")
  endif()
endif()

include(embed_resources)
include(GenerateExportHeader)

#####################################################
################# gathering info for config.h #######
#####################################################

feature_summary(WHAT ENABLED_FEATURES VAR feature_list)
string(REPLACE "\n" "\\n" feature_list "${feature_list}")

get_directory_property( compiler_definitions COMPILE_DEFINITIONS )
set(EXTRA_CXX_FLAGS_FROM_DEFS "")
foreach(DEF ${compiler_definitions})
  set(EXTRA_CXX_FLAGS_FROM_DEFS "${EXTRA_CXX_FLAGS_FROM_DEFS} -D${DEF}")
endforeach()

#####################################################
################# include the rest ##################
#####################################################
add_subdirectory(external_packages)
add_subdirectory(casadi)
add_subdirectory(experimental EXCLUDE_FROM_ALL)
add_subdirectory(misc)

option(WITH_EXAMPLES "Build examples" ON)
if(WITH_EXAMPLES)
  add_subdirectory(docs/examples)
  add_subdirectory(docs/api/examples/ctemplate)
endif()

#####################################################
######################### docs ######################
#####################################################
option(WITH_DOC "Enable documentation generation" OFF)
if(WITH_DOC)
  set(DOXYFILE_LATEX "NO")
  set(DOXYFILE_IN_DIR ${PROJECT_SOURCE_DIR}/docs/api)
  set(DOXYFILE_OUTPUT_DIR  ${PROJECT_SOURCE_DIR}/docs/api)
  set(DOXYFILE_SOURCE_DIR ${PROJECT_SOURCE_DIR}/casadi)
  include(UseDoxygen OPTIONAL)
  add_custom_target(tutor)
endif()

# Print a summary
feature_summary(WHAT ALL)

# what is this?
if(EXTRA_CMAKE)
  include(${EXTRA_CMAKE})
endif()

########################################################
########### Generate config files ######################
########################################################
# Add all targets to the build-tree export set
export(TARGETS casadi
  FILE "${PROJECT_BINARY_DIR}/casadi-targets.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE CasADi)

# Create casadi-config.cmake
file(RELATIVE_PATH REL_LIB_PREFIX "${CMAKE_PREFIX}" "${LIB_PREFIX}")
# ... for the build tree
set(CONF_LIBRARY_DIRS "${LIBRARY_OUTPUT_PATH}")
configure_file(casadi-config.cmake.in "${PROJECT_BINARY_DIR}/casadi-config.cmake" @ONLY)
# ... for the install tree
set(CONF_LIBRARY_DIRS "\${CASADI_CMAKE_DIR}/${REL_LIB_PREFIX}")
configure_file(casadi-config.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/casadi-config.cmake" @ONLY)
# Create casadi-config-version.cmake, same for both build and install tree
configure_file(casadi-config-version.cmake.in
  "${PROJECT_BINARY_DIR}/casadi-config-version.cmake" @ONLY)

# Install the casadi-config.cmake and casadi-config-version.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/casadi-config.cmake"
  "${PROJECT_BINARY_DIR}/casadi-config-version.cmake"
  DESTINATION "${CMAKE_PREFIX}" COMPONENT dev)

# Generate the casadi.pc.in configuration file
configure_file(${PROJECT_SOURCE_DIR}/misc/casadi.pc.in
  "${PROJECT_BINARY_DIR}/casadi.pc" @ONLY)

# Install the casadi.pc
install(FILES
  "${PROJECT_BINARY_DIR}/casadi.pc"
  DESTINATION "${LIB_PREFIX}/pkgconfig" COMPONENT dev)

# Install the export set for use with the install-tree
file(RELATIVE_PATH REL_CMAKE_PREFIX "${CMAKE_INSTALL_PREFIX}" "${CMAKE_PREFIX}")
install(EXPORT casadi-targets DESTINATION ${REL_CMAKE_PREFIX})

set(CPACK_PACKAGE_CONTACT "casadi-users@googlegroups.com")
set(CPACK_PACKAGE_VERSION ${PACKAGE_VERSION_FULL})
include(CPack)
