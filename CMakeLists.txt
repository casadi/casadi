cmake_minimum_required(VERSION 3.10.2)
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

set(CASADI_TP_SHARED_PREFIX casadi-tp-)

############################################################
####################### Policies ###########################
############################################################
# Enforce cmake2.4+ style behaviour for
# mixed filename/full path linker directives
cmake_policy(SET CMP0003 NEW)

# cmake_policy propagation rules
if(POLICY CMP0011)
cmake_policy(SET CMP0011 NEW)
endif()

# cmake_policy avoid escaping in add_definitions
cmake_policy(SET CMP0005 NEW)

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
set(BUILD_BLASFEO_VERSION "edf92b396adddd9e548b9786f87ad290a0971329" CACHE STRING "Tag/branch/hash to be used for BUILD_BLASFEO_GIT_REPO")
set(BUILD_BLASFEO_GIT_REPO "https://github.com/giaf/blasfeo.git" CACHE STRING "Git repo used together with BUILD_BLASFEO_VERSION")
option(WITH_BUILD_BLASFEO "Build BLASFEO (BUILD_BLASFEO_VERSION=${BUILD_BLASFEO_VERSION}) from downloaded source (BUILD_BLASFEO_GIT_REPO=${BUILD_BLASFEO_GIT_REPO})." OFF)

# HPIPM
option(WITH_HPIPM "Compile the interface to HPIPM" OFF)
set(BUILD_HPIPM_VERSION "0e0c9f4e0d4081dceafa9b37c396db50bce0e81a" CACHE STRING "Tag/branch/hash to be used for BUILD_HPIPM_GIT_REPO")
set(BUILD_HPIPM_GIT_REPO "https://github.com/jgillis/hpipm.git" CACHE STRING "Git repo used together with BUILD_HPIPM_VERSION")
option(WITH_BUILD_HPIPM "Build HPIPM (BUILD_HPIPM_VERSION=${BUILD_HPIPM_VERSION}) from downloaded source (BUILD_HPIPM_GIT_REPO=${BUILD_HPIPM_GIT_REPO})." OFF)
set(BUILD_HPIPM_DEPENDENCIES BLASFEO)

# SuperSCS: An conic solver
option(WITH_SUPERSCS "Compile the interface to SuperSCS (the source code for SuperSCS is included)" ${WITH_SUPERSCS_DEF})
set(BUILD_SUPERSCS_VERSION "4d2d1bd03ed4cf93e684a880b233760ce34ca69c" CACHE STRING "Tag/branch/hash to be used for BUILD_SUPERSCS_GIT_REPO")
set(BUILD_SUPERSCS_GIT_REPO "https://github.com/jgillis/scs.git" CACHE STRING "Git repo used together with BUILD_SUPERSCS_VERSION")
option(WITH_BUILD_SUPERSCS "Build SuperSCS (BUILD_SUPERSCS_VERSION=${BUILD_SUPERSCS_VERSION}) from downloaded source (BUILD_SUPERSCS_GIT_REPO=${BUILD_SUPERSCS_GIT_REPO})." OFF)
add_feature_info(superscs-interface WITH_SUPERSCS "Interface to Conic solver SuperSCS.")

# OSQP: An operator-splitting QP solver
option(WITH_OSQP "Compile the interface to OSQP" OFF)
set(BUILD_OSQP_VERSION "v0.6.2" CACHE STRING "Tag/branch/hash to be used for BUILD_OSQP_GIT_REPO")
set(BUILD_OSQP_GIT_REPO "https://github.com/osqp/osqp.git" CACHE STRING "Git repo used together with BUILD_OSQP_VERSION")
option(WITH_BUILD_OSQP "Build OSQP (BUILD_OSQP_VERSION=${BUILD_OSQP_VERSION}) from downloaded source (BUILD_OSQP_GIT_REPO=${BUILD_OSQP_GIT_REPO})." OFF)
add_feature_info(osqp-interface WITH_OSQP "Interface to QP solver OSQP.")

set(BUILD_EIGEN3_VERSION "3.4.0" CACHE STRING "Tag/branch/hash to be used for BUILD_EIGEN3_GIT_REPO")
set(BUILD_EIGEN3_GIT_REPO "https://gitlab.com/libeigen/eigen.git" CACHE STRING "Git repo used together with BUILD_EIGEN3_VERSION")
option(WITH_BUILD_EIGEN3 "Build Eigen (BUILD_EIGEN3_VERSION=${BUILD_EIGEN3_VERSION}) from downloaded source (BUILD_EIGEN3_GIT_REPO=${BUILD_EIGEN3_GIT_REPO})." OFF)

set(BUILD_SIMDE_VERSION "v0.7.2" CACHE STRING "Tag/branch/hash to be used for BUILD_SIMDE_GIT_REPO")
set(BUILD_SIMDE_GIT_REPO "https://github.com/simd-everywhere/simde.git" CACHE STRING "Git repo used together with BUILD_BLASFEO_VERSION")
option(WITH_BUILD_SIMDE "Build Simde (BUILD_SIMDE_VERSION=${BUILD_SIMDE_VERSION}) from downloaded source (BUILD_SIMDE_GIT_REPO=${BUILD_SIMDE_GIT_REPO})." OFF)


# PROXQP: A QP solver based on proximal methods
option(WITH_PROXQP "Compile the interface to PROXQP" OFF)
set(BUILD_PROXQP_VERSION "v0.3.2" CACHE STRING "Tag/branch/hash to be used for BUILD_PROXQP_GIT_REPO")
set(BUILD_PROXQP_GIT_REPO "https://github.com/Simple-Robotics/proxsuite.git" CACHE STRING "Git repo used together with BUILD_PROXQP_VERSION")
option(WITH_BUILD_PROXQP "Build PROXQP (BUILD_PROXQP_VERSION=${BUILD_PROXQP_VERSION}) from downloaded source (BUILD_PROXQP_GIT_REPO=${BUILD_PROXQP_GIT_REPO})." OFF)
set(BUILD_PROXQP_DEPENDENCIES EIGEN3 SIMDE)

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
set(BUILD_LAPACK_GIT_REPO https://github.com/xianyi/OpenBLAS.git  CACHE STRING "Git repo used together with BUILD_LAPACK_VERSION")
set(BUILD_LAPACK_VERSION v0.3.21 CACHE STRING "Tag/branch/hash to be used for BUILD_LAPACK_GIT_REPO")

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
set(BUILD_IPOPT_GIT_REPO "https://github.com/jgillis/Ipopt-1.git" CACHE STRING "Git repo used together with BUILD_IPOPT_VERSION")
set(BUILD_IPOPT_VERSION "3.14.11.mod" CACHE STRING "Tag/branch/hash to be used for BUILD_IPOPT_GIT_REPO")
option(WITH_BUILD_IPOPT "Build IPOPT (BUILD_IPOPT_VERSION=${BUILD_IPOPT_VERSION}) from downloaded source (BUILD_IPOPT_GIT_REPO=${BUILD_IPOPT_GIT_REPO})." OFF)
set(BUILD_IPOPT_DEPENDENCIES MUMPS)


option(WITH_MOCKUP_REQUIRED ON)
set(BUILD_MOCKUPS_GIT_REPO https://github.com/casadi/mockups.git CACHE STRING "Git repo used together with BUILD_MOCKUPS_VERSION")
set(BUILD_MOCKUPS_VERSION v57 CACHE STRING "Tag/branch/hash to be used for BUILD_MOCKUPS_GIT_REPO")

# Knitro: An NLP solver
option(WITH_KNITRO "Compile the interface to KNITRO" OFF)
option(WITH_MOCKUP_KNITRO "Use mockup KNITRO (BUILD_MOCKUPS_VERSION=${BUILD_MOCKUPS_VERSION}) from downloaded source (BUILD_MOCKUPS_GIT_REPO=${BUILD_MOCKUPS_GIT_REPO})." OFF)


# Worhp: A sparse SQP method
option(WITH_WORHP "Compile the WORHP interface" OFF)
option(WITH_MOCKUP_WORHP "Use mockup WORHP (BUILD_MOCKUPS_VERSION=${BUILD_MOCKUPS_VERSION}) from downloaded source (BUILD_MOCKUPS_GIT_REPO=${BUILD_MOCKUPS_GIT_REPO})." OFF)

# CPLEX: A sparse mixed-integer QP solver
option(WITH_CPLEX "Compile the interface to CPLEX" OFF)
option(WITH_CPLEX_SHARED "Use shared cplex libraries" OFF)
option(WITH_MOCKUP_CPLEX "Use mockup CPLEX (BUILD_MOCKUPS_VERSION=${BUILD_MOCKUPS_VERSION}) from downloaded source (BUILD_MOCKUPS_GIT_REPO=${BUILD_MOCKUPS_GIT_REPO})." OFF)

# GUROBI: A sparse mixed-integer QP solver
option(WITH_GUROBI "Compile the interface to GUROBI" OFF)
option(WITH_MOCKUP_GUROBI "Use mockup GUROBI (BUILD_MOCKUPS_VERSION=${BUILD_MOCKUPS_VERSION}) from downloaded source (BUILD_MOCKUPS_GIT_REPO=${BUILD_MOCKUPS_GIT_REPO})." OFF)

# Bonmin: An MINLP solver
option(WITH_BONMIN "Compile the interface to BONMIN" OFF)
set(BUILD_BONMIN_GIT_REPO https://github.com/coin-or/Bonmin.git CACHE STRING "Git repo used together with BUILD_BONMIN_VERSION")
set(BUILD_BONMIN_VERSION "releases/1.8.8" CACHE STRING "Tag/branch/hash to be used for BUILD_BONMIN_GIT_REPO")
option(WITH_BUILD_BONMIN "Build BONMIN (BUILD_BONMIN_VERSION=${BUILD_BONMIN_VERSION}) from downloaded source (BUILD_BONMIN_GIT_REPO=${BUILD_BONMIN_GIT_REPO})." OFF)
set(BUILD_BONMIN_DEPENDENCIES IPOPT CBC)

# CBC: A sparse mixed-integer linear solver
option(WITH_CBC "Compile the CBC interface" OFF)
set(BUILD_CBC_VERSION "releases/2.10.6" CACHE STRING "Tag/branch/hash to be used for BUILD_CBC_GIT_REPO")
set(BUILD_CBC_GIT_REPO "https://github.com/coin-or/Cbc.git" CACHE STRING "Git repo used together with BUILD_CBC_VERSION")

option(WITH_BUILD_CBC "Build CBC (BUILD_CBC_VERSION=${BUILD_CBC_VERSION}) from downloaded source." OFF)
set(BUILD_CBC_DEPENDENCIES CLP)

# CLP: A sparse linear solver
option(WITH_CLP "Compile the CLP interface" OFF)
set(BUILD_CLP_VERSION "releases/1.17.7" CACHE STRING "Tag/branch/hash to be used for BUILD_CLP_GIT_REPO")
set(BUILD_CLP_GIT_REPO "https://github.com/coin-or/Clp.git" CACHE STRING "Git repo used together with BUILD_CLP_VERSION")
set(BUILD_COINUTILS_VERSION "releases/2.11.6" CACHE STRING "Tag/branch/hash to be used for BUILD_COINUTILS_GIT_REPO")
set(BUILD_COINUTILS_GIT_REPO "https://github.com/coin-or/CoinUtils.git" CACHE STRING "Git repo used together with BUILD_COINUTILS_VERSION")
set(BUILD_CGL_VERSION "releases/0.60.4" CACHE STRING "Tag/branch/hash to be used for BUILD_CGL_GIT_REPO")
set(BUILD_CGL_GIT_REPO "https://github.com/coin-or/Cgl.git" CACHE STRING "Git repo used together with BUILD_CGL_VERSION")
set(BUILD_OSI_VERSION "releases/0.108.7" CACHE STRING "Tag/branch/hash to be used for BUILD_OSI_GIT_REPO")
set(BUILD_OSI_GIT_REPO "https://github.com/coin-or/Osi.git" CACHE STRING "Git repo used together with BUILD_OSI_VERSION")

option(WITH_BUILD_CLP "Build CLP (BUILD_CLP_VERSION=${BUILD_CLP_VERSION}) from downloaded source (BUILD_CLP_GIT_REPO=${BUILD_CLP_GIT_REPO})." OFF)
set(BUILD_CLP_DEPENDENCIES MUMPS LAPACK)

# MUMPS: Sparse direct linear solver
option(WITH_MUMPS "Enable MUMPS interface" OFF)
set(BUILD_MUMPS_TP_VERSION "releases/3.0.2" CACHE STRING "Tag/branch/hash to be used for BUILD_MUMPS_TP_GIT_REPO")
set(BUILD_METIS_TP_VERSION "releases/2.0.0" CACHE STRING "Tag/branch/hash to be used for BUILD_BETIS_TP_GIT_REPO")
set(BUILD_METIS_TP_GIT_REPO "https://github.com/coin-or-tools/ThirdParty-Metis.git" CACHE STRING "Git repo used together with BUILD_METIS_TP_VERSION")
option(WITH_BUILD_METIS "Build Metis from downloaded source." OFF)
set(BUILD_MUMPS_TP_GIT_REPO "https://github.com/coin-or-tools/ThirdParty-Mumps.git" CACHE STRING "Git repo used together with BUILD_MUMPS_TP_VERSION")

option(WITH_BUILD_MUMPS "Build MUMPS (BUILD_MUMPS_TP_VERSION=${BUILD_MUMPS_TP_VERSION}) from downloaded source (BUILD_MUMPS_TP_GIT_REPO=${BUILD_MUMPS_TP_GIT_REPO})." OFF)
set(BUILD_MUMPS_DEPENDENCIES LAPACK METIS)

# SPRAL
option(WITH_SPRAL "Enable SPRAL interface" OFF)
set(BUILD_SPRAL_VERSION "d385d2c9e858366d257cafaaf05760ffa6543e26" CACHE STRING "Tag/branch/hash to be used for BUILD_SPRAL_GIT_REPO")
set(BUILD_SPRAL_GIT_REPO "https://github.com/ralna/spral.git" CACHE STRING "Git repo used together with BUILD_SPRAL_VERSION")
option(WITH_BUILD_SPRAL "Build SPRAL (BUILD_SPRAL_VERSION=${BUILD_SPRAL_VERSION}) from downloaded source (BUILD_SPRAL_GIT_REPO=${BUILD_SPRAL_GIT_REPO})." OFF)
set(BUILD_SPRAL_DEPENDENCIES LAPACK METIS)

# HSL: Sparse direct linear solvers
option(WITH_HSL "Enable HSL interface" OFF)
option(WITH_MOCKUP_HSL "Use mockup WORHP (BUILD_MOCKUPS_VERSION=${BUILD_MOCKUPS_VERSION}) from downloaded source (BUILD_MOCKUPS_GIT_REPO=${BUILD_MOCKUPS_GIT_REPO})." OFF)
option(WITH_BUILD_HSL "Build HSL from supplied source directory (BUILD_HSL_SRC=${BUILD_HSL_SRC})." OFF)
if(WITH_BUILD_HSL AND WITH_MOCKUP_HSL)
    message(FATAL_ERROR "WITH_BUILD_HSL and WITH_MOCKUP_HSL are mutually exclusive.")
endif()
set(BUILD_HSL_SRC "." CACHE STRING "Source loction of coinhsl")
set(BUILD_HSL_DEPENDENCIES LAPACK METIS)

# HIGHS: A MILP / QP solver
option(WITH_HIGHS "Compile the HiGHS interface" OFF)
set(BUILD_HIGHS_GIT_REPO "https://github.com/ERGO-Code/HiGHS" CACHE STRING "Git repo used together with BUILD_HIGHS_VERSION")
set(BUILD_HIGHS_VERSION "v1.4.1" CACHE STRING "Tag/branch/hash to be used for BUILD_HIGHS_GIT_REPO")
option(WITH_BUILD_HIGHS "Build HiGHS (BUILD_HIGHS_VERSION=${BUILD_HIGHS_VERSION}) from downloaded source (BUILD_HIGHS_GIT_REPO=${BUILD_HIGHS_GIT_REPO})." OFF)


option(ALLOW_DOCKER "System has docker installed" OFF)

function(libraries2flags libraries out)
# libtool cannot handle absoulte paths (i.e. it omits the library from the list of dependencies)
  set(local "")
  foreach(LIB ${libraries})
    if (LIB MATCHES "^-")
      # Pass -l, -L as-is
      set(local "${local} ${LIB}")
    else()
      if (LIB MATCHES "System/Library/Frameworks/Accelerate.framework")
        #set(${out} "-framework Accelerate" PARENT_SCOPE)
        set(local "${local} -framework Accelerate")
      else()
          get_filename_component(LIB_DIR "${LIB}" DIRECTORY)
          get_filename_component(LIB_NAME "${LIB}" NAME)
          string(REGEX REPLACE "${CMAKE_SHARED_LIBRARY_PREFIX}([A-Za-z0-9_-]+)${CMAKE_SHARED_LIBRARY_SUFFIX}" "\\1" LIB_NAME_MOD ${LIB_NAME})
          if(LIB_NAME STREQUAL LIB_NAME_MOD)
            set(local "${local} ${LIB}")
          else()
            set(local "${local} -L${LIB_DIR} -l${LIB_NAME_MOD}")
          endif()
      endif()
    endif()
    set(${out} "${local}" PARENT_SCOPE)
  endforeach()
endfunction()


macro(eager_update DEP_NAME)
    set(VERSION "VERSION:${BUILD_${DEP_NAME}_VERSION}\nREPO:${BUILD_${DEP_NAME}_GIT_REPO}")
    set(VERSION_FILE "${CMAKE_BINARY_DIR}/external_projects/versions/${DEP_NAME}_version.txt")
    # Get lowercase version of DEP_NAME
    string(TOLOWER "${DEP_NAME}" DEP_NAME_LOWERCASE)

    if(NOT EXISTS ${VERSION_FILE})
      file(WRITE ${VERSION_FILE} "${VERSION}")
    endif()


    # Read the version file
    file(READ ${VERSION_FILE} VERSION_FILE_CONTENTS)
    if(NOT VERSION STREQUAL VERSION_FILE_CONTENTS)
      message("${DEP_NAME} is out of date. Removing old version.")
      file(WRITE ${VERSION_FILE} "${VERSION}")
      file(GLOB TO_REMOVE LIST_DIRECTORIES true "${CMAKE_BINARY_DIR}/external_projects/src/${DEP_NAME_LOWERCASE}-*")
      if(TO_REMOVE)
       file(REMOVE_RECURSE ${TO_REMOVE})
      endif()
    endif()
endmacro()

set(WITH_BUILD_MOCKUPS OFF)
# Packages must be sorted from high to low in the dependency graph
foreach(PKG CPLEX KNITRO GUROBI WORHP HSL)
  string(TOLOWER ${PKG} PKG_LOWER)
  if(WITH_${PKG} AND NOT WITH_BUILD_${PKG})
    if(NOT WITH_MOCKUP_${PKG})
      get_property(PACKAGES_NOT_FOUND_BACKUP GLOBAL PROPERTY PACKAGES_NOT_FOUND)
      get_property(PACKAGES_FOUND_BACKUP GLOBAL PROPERTY PACKAGES_FOUND)
      find_package(${PKG})
    endif()
    if (WITH_MOCKUP_${PKG} OR (WITH_MOCKUP_REQUIRED AND NOT ${PKG}_FOUND))

      # Add it to FeatureSummary
      get_property(PKG_HELP CACHE WITH_MOCKUP_${PKG} PROPERTY HELPSTRING)
      set(WITH_MOCKUP_${PKG} ON CACHE BOOL "${PKG_HELP}" FORCE)
      add_feature_info(${PKG_LOWER}-mockup-build ON ${PKG_HELP})

      # Remove from PACKAGES_NOT_FOUND
      set_property(GLOBAL PROPERTY PACKAGES_NOT_FOUND ${PACKAGES_NOT_FOUND_BACKUP})
      set_property(GLOBAL PROPERTY PACKAGES_FOUND ${PACKAGES_FOUND_BACKUP})

      # We are going to need build dependencies
      foreach(DEP ${BUILD_${PKG}_DEPENDENCIES})
        set(WITH_${DEP} ON)
      endforeach()
      set(${PKG}_FOUND ON)
      set(WITH_BUILD_MOCKUPS ON)
      eager_update(MOCKUPS)
    endif()
    if(NOT ${PKG}_FOUND)
      message(FATAL_ERROR "${PKG} not found on your system yet required. Set WITH_MOCKUP_${PKG} or WITH_MOCKUP_REQUIRED ON to build a mockup version such that the CasADi build can proceed.")
    endif()
  endif()
endforeach()

# Packages must be sorted from high to low in the dependency graph
foreach(PKG HIGHS PROXQP OSQP SUPERSCS BONMIN IPOPT CBC CLP MUMPS SPRAL HSL METIS HPIPM BLASFEO LAPACK EIGEN3 SIMDE)
  string(TOLOWER ${PKG} PKG_LOWER)
  if(WITH_${PKG} AND NOT WITH_MOCKUP_${PKG})
    if(NOT WITH_BUILD_${PKG})
      get_property(PACKAGES_NOT_FOUND_BACKUP GLOBAL PROPERTY PACKAGES_NOT_FOUND)
      get_property(PACKAGES_FOUND_BACKUP GLOBAL PROPERTY PACKAGES_FOUND)
      find_package(${PKG})
    endif()
    if (WITH_BUILD_${PKG} OR (WITH_BUILD_REQUIRED AND NOT ${PKG}_FOUND))

      # Add it to FeatureSummary
      message("get_property(PKG_HELP CACHE WITH_BUILD_${PKG} PROPERTY HELPSTRING)")
      get_property(PKG_HELP CACHE WITH_BUILD_${PKG} PROPERTY HELPSTRING)
      set(WITH_BUILD_${PKG} ON CACHE BOOL "${PKG_HELP}" FORCE)
      add_feature_info(${PKG_LOWER}-sourcebuild ON ${PKG_HELP})

      # Remove from PACKAGES_NOT_FOUND
      set_property(GLOBAL PROPERTY PACKAGES_NOT_FOUND ${PACKAGES_NOT_FOUND_BACKUP})
      set_property(GLOBAL PROPERTY PACKAGES_FOUND ${PACKAGES_FOUND_BACKUP})

      # We are going to need build dependencies
      foreach(DEP ${BUILD_${PKG}_DEPENDENCIES})
        set(WITH_${DEP} ON)
      endforeach()
      set(${PKG}_FOUND ON)
      eager_update(${PKG})
    endif()
    if(NOT ${PKG}_FOUND)
      message(FATAL_ERROR "${PKG} not found on your system yet required. Set WITH_BUILD_${PKG} or WITH_BUILD_REQUIRED ON to build this required package from source.")
    endif()
  endif()
endforeach()

set(SHARED_LIBRARY_RELDIR lib)
set(SHARED_LIBRARY64_RELDIR lib64)
set(CMAKE_SHARED_LIBRARY_PREFIX_SHORT ${CMAKE_SHARED_LIBRARY_PREFIX})
if (WIN32)
    set(SHARED_LIBRARY_RELDIR bin)
    set(CMAKE_SHARED_LIBRARY_PREFIX_SHORT "")
else()
    if (CMAKE_SIZEOF_VOID_P EQUAL 4)
        set(SHARED_LIBRARY64_RELDIR lib)
    endif()
endif()



if (WITH_BUILD_MOCKUPS)
ExternalProject_Add(mockups-external
  GIT_REPOSITORY ${BUILD_MOCKUPS_GIT_REPO}
  GIT_TAG ${BUILD_MOCKUPS_VERSION}
  PREFIX "${CMAKE_BINARY_DIR}/external_projects"
  UPDATE_COMMAND ""
  BUILD_ALWAYS ON
  CMAKE_COMMAND "${EXTERNAL_CMAKE_COMMAND}"
  CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>/mockups)
endif()

if(WITH_CPLEX AND WITH_MOCKUP_CPLEX)
  add_library(cplex::cplex SHARED IMPORTED)
  add_dependencies(cplex::cplex mockups-external)
  set(CPLEX_ADAPTOR_LIB "${CMAKE_BINARY_DIR}/external_projects/mockups/cplex/${SHARED_LIBRARY_RELDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}cplex_adaptor${CMAKE_SHARED_LIBRARY_SUFFIX}")
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/mockups/cplex/include")
  set_target_properties(cplex::cplex PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/mockups/cplex/include"
    IMPORTED_LOCATION "${CPLEX_ADAPTOR_LIB}"
    IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/mockups/cplex/lib/${CMAKE_IMPORT_LIBRARY_PREFIX}cplex_adaptor${CMAKE_IMPORT_LIBRARY_SUFFIX}"
  )
  install(FILES
  ${CPLEX_ADAPTOR_LIB}
  DESTINATION "${LIB_PREFIX}")
endif()
add_feature_info(cplex-interface WITH_CPLEX "Interface to the QP solver CPLEX.")

if(WITH_GUROBI AND WITH_MOCKUP_GUROBI)
  add_library(gurobi::gurobi SHARED IMPORTED)
  add_dependencies(gurobi::gurobi mockups-external)
  set(GUROBI_ADAPTOR_LIB "${CMAKE_BINARY_DIR}/external_projects/mockups/gurobi/${SHARED_LIBRARY_RELDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}gurobi_adaptor${CMAKE_SHARED_LIBRARY_SUFFIX}")
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/mockups/gurobi/include")
  set_target_properties(gurobi::gurobi PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/mockups/gurobi/include"
    IMPORTED_LOCATION "${GUROBI_ADAPTOR_LIB}"
    IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/mockups/gurobi/lib/${CMAKE_IMPORT_LIBRARY_PREFIX}gurobi_adaptor${CMAKE_IMPORT_LIBRARY_SUFFIX}"
  )
  install(FILES
  ${GUROBI_ADAPTOR_LIB}
  DESTINATION "${LIB_PREFIX}")
endif()
add_feature_info(gurobi-interface WITH_GUROBI "Interface to the (mixed-integer) QP solver GUROBI")

if(WITH_KNITRO AND WITH_MOCKUP_KNITRO)
  add_library(knitro::knitro SHARED IMPORTED)
  add_dependencies(knitro::knitro mockups-external)
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/mockups/knitro/include")
  set_target_properties(knitro::knitro PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/mockups/knitro/include"
    IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/mockups/knitro/${SHARED_LIBRARY_RELDIR}/${CMAKE_SHARED_LIBRARY_PREFIX_SHORT}knitro${CMAKE_SHARED_LIBRARY_SUFFIX}"
    IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/mockups/knitro/lib/${CMAKE_SHARED_LIBRARY_PREFIX_SHORT}knitro${CMAKE_IMPORT_LIBRARY_SUFFIX}"
  )
endif()
add_feature_info(knitro-interface WITH_KNITRO "Interface to the NLP solver KNITRO.")

if(WITH_WORHP AND WITH_MOCKUP_WORHP)
  add_library(worhp::worhp SHARED IMPORTED)
  add_dependencies(worhp::worhp mockups-external)
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/mockups/worhp/include")
  set_target_properties(worhp::worhp PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/mockups/worhp/include"
    IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/mockups/worhp/${SHARED_LIBRARY_RELDIR}/${CMAKE_SHARED_LIBRARY_PREFIX_SHORT}worhp${CMAKE_SHARED_LIBRARY_SUFFIX}"
    IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/mockups/worhp/lib/${CMAKE_SHARED_LIBRARY_PREFIX_SHORT}worhp${CMAKE_IMPORT_LIBRARY_SUFFIX}"
  )
endif()
add_feature_info(worhp-interface WITH_WORHP "Interface to the NLP solver Worhp (requires gfortran, gomp).")

if(WITH_LAPACK)
  if (WITH_BUILD_LAPACK)
    # Requires gcc>=5.0
    ExternalProject_Add(openblas-external
      GIT_REPOSITORY ${BUILD_LAPACK_GIT_REPO}
      GIT_TAG ${BUILD_LAPACK_VERSION}
      PREFIX "${CMAKE_BINARY_DIR}/external_projects"
      UPDATE_COMMAND ""
      CMAKE_COMMAND "${EXTERNAL_CMAKE_COMMAND}"
      PATCH_COMMAND ${EXTERNAL_CMAKE_COMMAND} -DFILE=CMakeLists.txt -DCASADI_TP_SHARED_PREFIX=${CASADI_TP_SHARED_PREFIX} -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/patch_openblas.cmake
      CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_LIBDIR=lib -DTARGET=CORE2 -DUSE_THREAD=ON -DNUM_THREADS=16 -DBUILD_SINGLE=ON -DBUILD_DOUBLE=ON -DBUILD_COMPLEX=ON -DBUILD_COMPLEX16=ON -DDYNAMIC_ARCH=OFF -DUSE_OPENMP=OFF -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
    set(LAPACK_LIBRARIES "${CMAKE_BINARY_DIR}/external_projects/${SHARED_LIBRARY_RELDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}${CASADI_TP_SHARED_PREFIX}openblas${CMAKE_SHARED_LIBRARY_SUFFIX}")    
    set(LAPACK_DEPENDS_TARGET openblas-external)

    add_library(lapack SHARED IMPORTED)
    add_dependencies(lapack openblas-external)
    set_target_properties(lapack PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/${SHARED_LIBRARY_RELDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}${CASADI_TP_SHARED_PREFIX}openblas${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_IMPORT_LIBRARY_PREFIX}${CASADI_TP_SHARED_PREFIX}openblas${CMAKE_IMPORT_LIBRARY_SUFFIX}"
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
    set(C_FLAGS "")
    if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
        set(C_FLAGS "-DCMAKE_C_FLAGS=-std=c99")
    endif()
    ExternalProject_Add(blasfeo-external
      GIT_REPOSITORY ${BUILD_BLASFEO_GIT_REPO}
      GIT_TAG ${BUILD_BLASFEO_VERSION}
      PREFIX "${CMAKE_BINARY_DIR}/external_projects"
      UPDATE_COMMAND ""
      CMAKE_COMMAND "${EXTERNAL_CMAKE_COMMAND}"
      CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} ${C_FLAGS} -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
    set(BLASFEO_PATH "${CMAKE_BINARY_DIR}/external_projects")
    set(BLASFEO_DEPENDS_TARGET blasfeo-external)

    add_library(blasfeo SHARED IMPORTED)
    add_dependencies(blasfeo blasfeo-external)
    set_target_properties(blasfeo PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/include"
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/${SHARED_LIBRARY_RELDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}blasfeo${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_IMPORT_LIBRARY_PREFIX}blasfeo${CMAKE_IMPORT_LIBRARY_SUFFIX}"
    )
    file(GLOB BLASFEO_HEADERS "${CMAKE_BINARY_DIR}/external_projects/include/blasfeo*")
    install(FILES
      ${BLASFEO_HEADERS}
      DESTINATION "${INCLUDE_PREFIX}")
  else()
    #add_library(lapack INTERFACE)
    #target_link_libraries(lapack INTERFACE ${LAPACK_LIBRARIES})
  endif()
endif()

message("CMAKE_C_COMPILER_ID: ${CMAKE_C_COMPILER_ID}")
if(WITH_HPIPM)
  if (WITH_BUILD_HPIPM)
    ExternalProject_Add(hpipm-external
      DEPENDS ${BLASFEO_DEPENDS_TARGET}
      GIT_REPOSITORY ${BUILD_HPIPM_GIT_REPO}
      GIT_TAG ${BUILD_HPIPM_VERSION}
      PREFIX "${CMAKE_BINARY_DIR}/external_projects"
      UPDATE_COMMAND ""
      CMAKE_COMMAND "${EXTERNAL_CMAKE_COMMAND}"
      CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DHPIPM_TESTING=OFF -DBUILD_SHARED_LIBS=ON -DBLASFEO_PATH=${BLASFEO_PATH} -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
    add_library(hpipm SHARED IMPORTED)
    add_dependencies(hpipm hpipm-external)
    target_link_libraries(hpipm INTERFACE blasfeo)
    set_target_properties(hpipm PROPERTIES
       IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/${SHARED_LIBRARY_RELDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}hpipm${CMAKE_SHARED_LIBRARY_SUFFIX}"
       IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_IMPORT_LIBRARY_PREFIX}hpipm${CMAKE_IMPORT_LIBRARY_SUFFIX}"
    )
    file(GLOB HPIPM_HEADERS "${CMAKE_BINARY_DIR}/external_projects/include/hpipm*")
    install(FILES
      ${HPIPM_HEADERS}
      DESTINATION "${INCLUDE_PREFIX}")
  else()
    #add_library(lapack INTERFACE)
    #target_link_libraries(lapack INTERFACE ${LAPACK_LIBRARIES})
  endif()
endif()

message(STATUS "CMAKE_SYSTEM_NAME: ${CMAKE_SYSTEM_NAME} ${WIN32}")

add_feature_info(lapack-interface WITH_LAPACK "Interface to LAPACK.")

if(WITH_BUILD_METIS)
  add_feature_info(metis-sourcebuild ON "Build METIS (BUILD_METIS_TP_VERSION=${BUILD_METIS_TP_VERSION}) from downloaded source.")
  if(ALLOW_DOCKER)
     set(COIN_RECONFIG_CMD docker run --rm -v <SOURCE_DIR>:/work ghcr.io/casadi/coinbuild:latest /bin/bash ./run_autotools /work COMMAND)
  else()
     set(COIN_RECONFIG_CMD)
  endif()
  ExternalProject_Add(metis-external
        GIT_REPOSITORY ${BUILD_METIS_TP_GIT_REPO}
        GIT_TAG ${BUILD_METIS_TP_VERSION}
        DEPENDS ${LAPACK_DEPENDS_TARGET}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND ${COIN_RECONFIG_CMD} <SOURCE_DIR>/configure ${CROSS_FLAGS} --prefix=<INSTALL_DIR>

          # GKFree is not defined (variadic definition in utils.c); fails on MacOS target
          "CFLAGS=-Wno-implicit-function-declaration"
        PATCH_COMMAND COMMAND cd <SOURCE_DIR> && ./get.Metis)
  add_library(metis INTERFACE IMPORTED)
  add_dependencies(metis metis-external)
  set_target_properties(metis PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/include/coin-or/metis"
    IMPORTED_LIBNAME "coinmetis"
  )
  message("metis exported")
endif()

if(WITH_HSL)
  if (WITH_MOCKUP_HSL)
      add_library(hsl::hsl SHARED IMPORTED)
      add_dependencies(hsl::hsl mockups-external)
      set_target_properties(hsl::hsl PROPERTIES
        IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/mockups/hsl/${SHARED_LIBRARY_RELDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}hsl${CMAKE_SHARED_LIBRARY_SUFFIX}"
        IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/mockups/hsl/lib/${CMAKE_IMPORT_LIBRARY_PREFIX}hsl${CMAKE_IMPORT_LIBRARY_SUFFIX}"
      )
  elseif(WITH_BUILD_HSL)
    add_feature_info(hsl-sourcebuild ON "Build CoinHSL from supplied source.")
    get_target_property(METIS_INCLUDE_DIR metis INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(METIS_LIB metis IMPORTED_LIBNAME)
    message("metis: ${METIS_LIB} ${METIS_INCLUDE_DIR}")
    # Copy source

    ExternalProject_Add(hsl-external
        SOURCE_DIR ${BUILD_HSL_SRC}
        DEPENDS metis
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        PATCH_COMMAND ${EXTERNAL_CMAKE_COMMAND} -DSRC=<SOURCE_DIR>  -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/patch_hsl.cmake
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --disable-static "--with-blas=${LAPACK_LFLAGS}" ${CROSS_FLAGS} "FCFLAGS=${LAPACK_LFLAGS} -L${CMAKE_BINARY_DIR}/external_projects/lib -Wl,-rpath,${CMAKE_BINARY_DIR}/external_projects/lib -l${METIS_LIB} -lm" --prefix=<INSTALL_DIR>
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        COMMAND cp ${CMAKE_BINARY_DIR}/external_projects/${SHARED_LIBRARY_RELDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}coinhsl${CMAKE_SHARED_LIBRARY_SUFFIX} ${CMAKE_BINARY_DIR}/external_projects/${SHARED_LIBRARY_RELDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}hsl${CMAKE_SHARED_LIBRARY_SUFFIX})
      add_library(hsl::hsl SHARED IMPORTED)
      add_dependencies(hsl::hsl hsl-external coinmetis lapack)
      set_target_properties(hsl::hsl PROPERTIES
        IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/${SHARED_LIBRARY_RELDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}coinhsl${CMAKE_SHARED_LIBRARY_SUFFIX}"
        IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_IMPORT_LIBRARY_PREFIX}coinhsl${CMAKE_IMPORT_LIBRARY_SUFFIX}"
      )
  endif()
endif()

if(WITH_MUMPS)
  if(WITH_BUILD_MUMPS)
    get_target_property(METIS_INCLUDE_DIR metis INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(METIS_LIB metis IMPORTED_LIBNAME)
    ExternalProject_Add(mumps-external
        GIT_REPOSITORY ${BUILD_MUMPS_TP_GIT_REPO}
        GIT_TAG ${BUILD_MUMPS_TP_VERSION}
        DEPENDS ${LAPACK_DEPENDS_TARGET} metis
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS}  --prefix=<INSTALL_DIR> "--with-metis-cflags=-I${METIS_INCLUDE_DIR}" "--with-metis-lflags=-L${CMAKE_BINARY_DIR}/external_projects/lib -l${METIS_LIB} -lm" "--with-lapack-lflags=${LAPACK_LFLAGS}" ${CROSS_FLAGS_FORTRAN}
        PATCH_COMMAND COMMAND cd <SOURCE_DIR> && ./get.Mumps)

    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/include/coin-or")
    add_library(mumps SHARED IMPORTED)
    add_dependencies(mumps mumps-external)
    set_target_properties(mumps PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}coinmumps${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}coinmumps${CMAKE_IMPORT_LIBRARY_SUFFIX}"
      INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/include/coin-or"
    )
  endif()
endif()
add_feature_info(mumps-interface WITH_MUMPS "Interface to MUMPS.")

set(SPRAL_DEPENDS_TARGET )
set(SPRAL_FLAGS )

message("WITH_SPRAL ${WITH_SPRAL} WITH_BUILD_SPRAL ${WITH_BUILD_SPRAL}")
if(WITH_SPRAL)

  # configure performed to resolve
  #      PKG_PROG_PKG_CONFIG: command not found"
  #      syntax error near unexpected token `HWLOC,'
  #      `PKG_CHECK_MODULES(HWLOC, hwloc,'
  if(ALLOW_DOCKER)
    set(SPRAL_RECONFIG_CMD docker run --rm -v <SOURCE_DIR>:/work ghcr.io/casadi/coinbuild:latest autoreconf -ivf -I /usr/share/aclocal /work COMMAND)
  else()
    set(SPRAL_RECONFIG_CMD autoreconf -ivf -I /usr/share/aclocal <SOURCE_DIR> COMMAND)
  endif()
  if(WITH_BUILD_SPRAL)
    get_target_property(METIS_INCLUDE_DIR metis INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(METIS_LIB metis IMPORTED_LIBNAME)
    ExternalProject_Add(spral-external
        GIT_REPOSITORY ${BUILD_SPRAL_GIT_REPO}
        GIT_TAG ${BUILD_SPRAL_VERSION}
        DEPENDS ${LAPACK_DEPENDS_TARGET} metis
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND ${SPRAL_RECONFIG_CMD} CFLAGS=-fPIC CXXFLAGS=-fPIC FFLAGS=-fPIC FCFLAGS=-fPIC <SOURCE_DIR>/configure ${CROSS_FLAGS} --disable-gpu --prefix=<INSTALL_DIR> --disable-openmp "--with-metis-inc-dir=${METIS_INCLUDE_DIR}" "--with-metis=-L${CMAKE_BINARY_DIR}/external_projects/lib -l${METIS_LIB} -lm" "--with-lapack=${LAPACK_LFLAGS}" "--with-blas=${LAPACK_LFLAGS}" ${CROSS_FLAGS_FORTRAN})
    set(SPRAL_DEPENDS_TARGET spral-external)
    set(SPRAL_FLAGS "--with-spral-cflags=-I${CMAKE_BINARY_DIR}/external_projects/include"
        "--with-spral-lflags=-L${CMAKE_BINARY_DIR}/external_projects/lib -lspral -lgfortran -lm -lcoinmetis ${LAPACK_LFLAGS}")
    install(FILES ${CMAKE_BINARY_DIR}/external_projects/src/spral-external/LICENCE DESTINATION ${INCLUDE_PREFIX}/licenses RENAME spral.txt)
  endif()
endif()
add_feature_info(spral-interface WITH_SPRAL "Interface to SPRAL.")

if((WITH_CLP AND WITH_BUILD_CLP) OR (WITH_CBC AND WITH_BUILD_CBC))
    add_feature_info(coinutils-sourcebuild ON "Build COINUTILS (BUILD_COINUTILS_VERSION=${BUILD_COINUTILS_VERSION}) from downloaded source.")
    ExternalProject_Add(coinutils-external
          DEPENDS ${LAPACK_DEPENDS_TARGET}
          GIT_REPOSITORY ${BUILD_COINUTILS_GIT_REPO}
          GIT_TAG ${BUILD_COINUTILS_VERSION}
          PREFIX "${CMAKE_BINARY_DIR}/external_projects"
          UPDATE_COMMAND ""
          CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS} "--with-lapack-lflags=${LAPACK_LFLAGS}" --prefix=<INSTALL_DIR> )
    add_library(coinutils SHARED IMPORTED)
    add_dependencies(coinutils coinutils-external)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/include/coin")
    set_target_properties(coinutils PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}CoinUtils${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}CoinUtils${CMAKE_IMPORT_LIBRARY_SUFFIX}"
    )
endif()

if(WITH_CLP)
  if(WITH_BUILD_CLP)
    add_feature_info(osi-sourcebuild ON "Build OSI (BUILD_OSI_VERSION=${BUILD_OSI_VERSION}) from downloaded source.")
    ExternalProject_Add(osi-external
        DEPENDS coinutils-external
        GIT_REPOSITORY ${BUILD_OSI_GIT_REPO}
        GIT_TAG ${BUILD_OSI_VERSION}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS} --prefix=<INSTALL_DIR> ${CROSS_FLAGS_FORTRAN}
	# For unknown reasons, Coin-or/configure disables pkg-config by default for clang compilers 
	--enable-pkg-config
    )
    add_library(osi SHARED IMPORTED)
    add_dependencies(osi osi-external)
    set_target_properties(osi PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}Osi${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}Osi${CMAKE_IMPORT_LIBRARY_SUFFIX}"
    )
    set(APPLE_FLAGS "")
    if(APPLE)
        # Delocate step wants to grow zlib install name
        set(APPLE_FLAGS "ADD_CXXFLAGS=-Wl,-headerpad_max_install_names")
    endif()
    ExternalProject_Add(clp-external
        DEPENDS osi-external mumps-external coinutils-external ${LAPACK_DEPENDS_TARGET}
        GIT_REPOSITORY ${BUILD_CLP_GIT_REPO}
        GIT_TAG ${BUILD_CLP_VERSION}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS} "--with-lapack-lflags=${LAPACK_LFLAGS}" --prefix=<INSTALL_DIR>  ${CROSS_FLAGS_FORTRAN}
	# For unknown reasons, Coin-or/configure disables pkg-config by default for clang compilers 
	--enable-pkg-config
	${APPLE_FLAGS}
	)

    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/include/coin")
    add_library(clp-raw SHARED IMPORTED)
    add_dependencies(clp-raw clp-external)
    set_target_properties(clp-raw PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}Clp${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}Clp${CMAKE_IMPORT_LIBRARY_SUFFIX}"
      INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/include/coin"
    )
    add_library(osi-clp SHARED IMPORTED)
    add_dependencies(osi-clp clp-external)
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
        GIT_REPOSITORY ${BUILD_CGL_GIT_REPO}
        GIT_TAG ${BUILD_CGL_VERSION}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS} --prefix=<INSTALL_DIR> ${CROSS_FLAGS_FORTRAN}
	# For unknown reasons, Coin-or/configure disables pkg-config by default for clang compilers 
	--enable-pkg-config
	)

    ExternalProject_Add(cbc-external
        DEPENDS coinutils-external cgl-external osi-external mumps-external clp-external ${LAPACK_DEPENDS_TARGET}
        GIT_REPOSITORY ${BUILD_CBC_GIT_REPO}
        GIT_TAG ${BUILD_CBC_VERSION}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS} --prefix=<INSTALL_DIR> "--with-lapack-lflags=${LAPACK_LFLAGS}"  ${CROSS_FLAGS_FORTRAN}
	# For unknown reasons, Coin-or/configure disables pkg-config by default for clang compilers 
	--enable-pkg-config
	)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/include/coin")
    add_library(cbc-raw SHARED IMPORTED)
    add_dependencies(cbc-raw cbc-external)
    set_target_properties(cbc-raw PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}Cbc${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}Cbc${CMAKE_IMPORT_LIBRARY_SUFFIX}"
      INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/include/coin"
    )
    add_library(cbc-solver SHARED IMPORTED)
    add_dependencies(cbc-solver cbc-external)
    set_target_properties(cbc-solver PROPERTIES
      IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}CbcSolver${CMAKE_SHARED_LIBRARY_SUFFIX}"
      IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}CbcSolver${CMAKE_IMPORT_LIBRARY_SUFFIX}"
    )
    add_library(cbc INTERFACE)
    target_link_libraries(cbc INTERFACE cbc-raw cbc-solver coinutils osi osi-clp)
  endif()
endif()
add_feature_info(cbc-interface WITH_CBC "Interface to the LP solver CBC.")

set(COIN_FLAGS "")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
set(COIN_FLAGS "--enable-debug")
endif()


if(WITH_IPOPT)
  if(WITH_BUILD_IPOPT)
    ExternalProject_Add(ipopt-external
        DEPENDS ${LAPACK_DEPENDS_TARGET} mumps-external ${SPRAL_DEPENDS_TARGET}
        GIT_REPOSITORY ${BUILD_IPOPT_GIT_REPO}
        GIT_TAG ${BUILD_IPOPT_VERSION}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS} --prefix=<INSTALL_DIR> "--with-mumps-cflags=-I${CMAKE_BINARY_DIR}/external_projects/include/coin-or/mumps" "--with-mumps-lflags=-L${CMAKE_BINARY_DIR}/external_projects/lib ${WHOLE_ARCHIVE_ON}-lcoinmumps ${WHOLE_ARCHIVE_OFF}"
        ${SPRAL_FLAGS}
        "--with-lapack-lflags= ${LAPACK_LFLAGS}" ${CROSS_FLAGS_FORTRAN}
	--disable-java
	# For unknown reasons, Coin-or/configure disables pkg-config by default for clang compilers 
	--enable-pkg-config
	--enable-dependency-linking
	${COIN_FLAGS}
	)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/include/coin-or")
    set(WITH_IPOPT_CALLBACK ON)

    add_library(ipopt SHARED IMPORTED)
    add_dependencies(ipopt ipopt-external)
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
        GIT_REPOSITORY ${BUILD_BONMIN_GIT_REPO}
        GIT_TAG ${BUILD_BONMIN_VERSION}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${CROSS_FLAGS} --prefix=<INSTALL_DIR>  ${CROSS_FLAGS_FORTRAN}
	# For unknown reasons, Coin-or/configure disables pkg-config by default for clang compilers 
	--enable-pkg-config
	# issue #2988
	--enable-dependency-linking
	${MINGW_FLAG}
	)

    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/include/coin")
    add_library(bonmin-raw SHARED IMPORTED)
    add_dependencies(bonmin-raw bonmin-external)
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

if(WITH_HIGHS)
  if(WITH_BUILD_HIGHS)
    ExternalProject_Add(highs-external
        GIT_REPOSITORY ${BUILD_HIGHS_GIT_REPO}
        GIT_TAG ${BUILD_HIGHS_VERSION}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CMAKE_COMMAND "${EXTERNAL_CMAKE_COMMAND}"
        CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/include/highs")
    add_library(highs SHARED IMPORTED)
    add_dependencies(highs highs-external)
    set_target_properties(highs PROPERTIES
       INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/include/highs"
       IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}highs${CMAKE_SHARED_LIBRARY_SUFFIX}"
       IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_IMPORT_LIBRARY_PREFIX}highs${CMAKE_IMPORT_LIBRARY_SUFFIX}"
    )
    install(DIRECTORY
      "${CMAKE_BINARY_DIR}/external_projects/include/highs"
      DESTINATION "${INCLUDE_PREFIX}")
  endif()
endif()
add_feature_info(highs-interface WITH_HIGHS "Interface to the MILP / QP solver HiGHS.")

if(WITH_BUILD_EIGEN3)
  ExternalProject_Add(eigen3-external
      GIT_REPOSITORY ${BUILD_EIGEN3_GIT_REPO}
      GIT_TAG ${BUILD_EIGEN3_VERSION}
      PREFIX "${CMAKE_BINARY_DIR}/external_projects"
      UPDATE_COMMAND ""
      CMAKE_COMMAND "${EXTERNAL_CMAKE_COMMAND}"
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
  add_library(eigen3 INTERFACE)
  add_dependencies(eigen3 eigen3-external)
  target_include_directories(eigen3 INTERFACE "${CMAKE_BINARY_DIR}/external_projects/include/eigen3")
  install(DIRECTORY
    "${CMAKE_BINARY_DIR}/external_projects/include/eigen3"
    DESTINATION "${INCLUDE_PREFIX}")
endif()

if(WITH_BUILD_SIMDE)
  ExternalProject_Add(simde-external
      GIT_REPOSITORY ${BUILD_SIMDE_GIT_REPO}
      GIT_TAG ${BUILD_SIMDE_VERSION}
      PREFIX "${CMAKE_BINARY_DIR}/external_projects"
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ${EXTERNAL_CMAKE_COMMAND} -E copy_directory
                    <SOURCE_DIR>/simde <INSTALL_DIR>/include/simde)
  add_library(Simde INTERFACE)
  add_dependencies(Simde Simde-external)
  target_include_directories(Simde INTERFACE "${CMAKE_BINARY_DIR}/external_projects/include/simde")
  install(DIRECTORY
    "${CMAKE_BINARY_DIR}/external_projects/include/simde"
    DESTINATION "${INCLUDE_PREFIX}")
endif()

if(WITH_PROXQP)
  if(WITH_BUILD_PROXQP)
    ExternalProject_Add(proxqp-external
        DEPENDS Simde eigen3
        GIT_REPOSITORY ${BUILD_PROXQP_GIT_REPO}
        GIT_TAG ${BUILD_PROXQP_VERSION}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CMAKE_COMMAND "${EXTERNAL_CMAKE_COMMAND}"
        CMAKE_ARGS --debug-find -DCMAKE_FIND_ROOT_PATH_MODE_PACKAGE=BOTH -DCMAKE_FIND_ROOT_PATH_MODE_INCLUDE=BOTH -DSimde_DIR=${CMAKE_BINARY_DIR}/external_projects/share/eigen3/cmake -DEigen3_DIR=${CMAKE_BINARY_DIR}/external_projects/share/eigen3/cmake -DCMAKE_PREFIX_PATH=${CMAKE_BINARY_DIR}/external_projects/include/ -DCMAKE_INCLUDE_PATH=${CMAKE_BINARY_DIR}/external_projects/include/ -DBUILD_TESTING=OFF -DBUILD_PYTHON_INTERFACE=OFF -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
        # --debug-find
    add_library(proxqp INTERFACE)
    add_dependencies(proxqp proxqp-external)
    target_include_directories(proxqp INTERFACE "${CMAKE_BINARY_DIR}/external_projects/include")
    target_link_libraries(proxqp INTERFACE eigen3)
    install(DIRECTORY
      "${CMAKE_BINARY_DIR}/external_projects/include/proxsuite"
      DESTINATION "${INCLUDE_PREFIX}")
  endif()
endif()
add_feature_info(proxqp-interface WITH_PROXQP "Interface to QP solver PROXQP.")

if(WITH_OSQP)
  if(WITH_BUILD_OSQP)
    ExternalProject_Add(osqp-external
        GIT_REPOSITORY ${BUILD_OSQP_GIT_REPO}
        GIT_TAG ${BUILD_OSQP_VERSION}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CMAKE_COMMAND "${EXTERNAL_CMAKE_COMMAND}"
        CMAKE_ARGS -DCMAKE_C_STANDARD_REQUIRED=TRUE -DCMAKE_C_STANDARD=99 -DENABLE_MKL_PARDISO=OFF -DPROFILING=OFF -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/include/osqp")
    add_library(osqp SHARED IMPORTED)
    add_dependencies(osqp osqp-external)
    set_target_properties(osqp PROPERTIES
       INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/include/osqp"
       IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}osqp${CMAKE_SHARED_LIBRARY_SUFFIX}"
       IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_IMPORT_LIBRARY_PREFIX}osqp${CMAKE_IMPORT_LIBRARY_SUFFIX}"
    )
    install(DIRECTORY
      "${CMAKE_BINARY_DIR}/external_projects/include/osqp"
      DESTINATION "${INCLUDE_PREFIX}")
    if(CMAKE_CROSS_COMPILING AND MINGW)
    add_custom_command(TARGET osqp-external POST_BUILD
      COMMAND gendef "${CMAKE_BINARY_DIR}/external_projects/${SHARED_LIBRARY_RELDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}osqp${CMAKE_SHARED_LIBRARY_SUFFIX}"
      COMMAND $ENV{DLLTOOL} --dllname "${CMAKE_BINARY_DIR}/external_projects/${SHARED_LIBRARY_RELDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}osqp${CMAKE_SHARED_LIBRARY_SUFFIX}" --def libosqp.def --output-lib osqp.lib
    )
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/osqp.lib"
      DESTINATION ${DEF_LIB_PREFIX}
    )
    endif()
  endif()
endif()
add_feature_info(osqp-interface WITH_OSQP "Interface to QP solver OSQP.")

if(WITH_SUPERSCS)
  if(WITH_BUILD_SUPERSCS)
    ExternalProject_Add(superscs-external
        GIT_REPOSITORY ${BUILD_SUPERSCS_GIT_REPO}
        GIT_TAG ${BUILD_SUPERSCS_VERSION}
        DEPENDS ${LAPACK_DEPENDS_TARGET}
        PREFIX "${CMAKE_BINARY_DIR}/external_projects"
        UPDATE_COMMAND ""
        CMAKE_COMMAND "${EXTERNAL_CMAKE_COMMAND}"
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external_projects/include/superscs")
    add_library(superscs STATIC IMPORTED)
    add_library(superscs-linsys STATIC IMPORTED)
    add_library(superscs-linsys-indirect STATIC IMPORTED)    
    add_dependencies(superscs superscs-external superscs-linsys superscs-linsys-indirect casadi)
    target_link_libraries(superscs INTERFACE lapack superscs-linsys superscs-linsys-indirect casadi)
    target_link_libraries(superscs-linsys-indirect INTERFACE superscs)
    set_target_properties(superscs PROPERTIES
       INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/include/superscs"
       IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}superscs${CMAKE_STATIC_LIBRARY_SUFFIX}"
    )
    set_target_properties(superscs-linsys PROPERTIES
       IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}linsys${CMAKE_STATIC_LIBRARY_SUFFIX}"
    )
    set_target_properties(superscs-linsys-indirect PROPERTIES
       IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_SHARED_LIBRARY_PREFIX}indirect${CMAKE_STATIC_LIBRARY_SUFFIX}"
    )
    install(DIRECTORY
      "${CMAKE_BINARY_DIR}/external_projects/include/superscs"
      DESTINATION "${INCLUDE_PREFIX}")
  endif()
endif()
add_feature_info(superscs-interface WITH_SUPERSCS "Interface to QP solver SUPERSCS.")

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

# OOQP: A sparse interior point QP solver
option(WITH_OOQP "Enable OOQP interface" OFF)
if(WITH_OOQP)
  find_package(OOQP REQUIRED)
endif()
add_feature_info(ooqp-interface WITH_OOQP "Interface to the QP solver OOQP (requires BLAS and HSL libraries).")

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
