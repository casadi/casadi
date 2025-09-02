cmake_minimum_required(VERSION 3.25...4.1)
set(CMAKE_CXX_SCAN_FOR_MODULES Off)
project(alpaqa
    VERSION 1.1.0
    DESCRIPTION "alpaqa augmented Lagrangian nonlinear programming solvers."
    HOMEPAGE_URL "https://github.com/kul-optec/alpaqa"
    LANGUAGES CXX
)
set(PY_VERSION_SUFFIX "a1")
include(CTest)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/find")

# Detect compiler features
include(cmake/DetectFeatures.cmake)

# Options
include(CMakeDependentOption)

# Enable/disable interfaces to other programming languages
option(ALPAQA_WITH_PYTHON
    "Build the Python bindings" Off)
option(ALPAQA_WITH_MATLAB
    "Build the MATLAB MEX interface" Off)
# Compiler support
option(ALPAQA_WITH_CXX_23
    "Enable components that require C++23" ${ALPAQA_WITH_CXX_23_DEFAULT})
# Enable/disable optional dependencies
option(ALPAQA_WITH_CASADI
    "Build the CasADi loader" On)
cmake_dependent_option(ALPAQA_WITH_EXTERNAL_CASADI
    "Link to the CasADi C++ library. When enabled, CasADi functions can be used to initialize alpaqa problems, otherwise, only pre-compiled CasADi functions can be loaded" Off "ALPAQA_WITH_CASADI" Off)
cmake_dependent_option(ALPAQA_WITH_CUTEST
    "Build the CUTEst loader" On "ALPAQA_WITH_CXX_23" Off)
option(ALPAQA_WITH_QPALM
    "Build the QPALM adapter" Off)
option(ALPAQA_WITH_IPOPT
    "Build the Ipopt adapter" Off)
option(ALPAQA_WITH_JSON
    "Enable JSON (de)serialization of options" Off)
cmake_dependent_option(ALPAQA_WITH_LBFGSB
    "Include the L-BFGS-B solver" On "ALPAQA_HAVE_FORTRAN" Off)
# Enable/disable optional components
option(ALPAQA_WITH_TESTS
    "Build the tests" ${BUILD_TESTING})
option(ALPAQA_WITH_EXAMPLES
    "Build the examples" On)
option(ALPAQA_WITH_PYTHON_PROBLEM_LOADER
    "Build the dynamic loader for problems defined in Python" Off)
option(ALPAQA_WITH_DRIVERS
    "Build the solver driver programs" On)
option(ALPAQA_WITH_GRADIENT_CHECKER
    "Build the utility for checking problem gradients" Off)
option(ALPAQA_WITH_OCP
    "Enable solvers tailored for optimal control problems" On)
cmake_dependent_option(ALPAQA_WITH_CASADI_OCP
    "Build the CasADi loader for optimal control problems"
    ${ALPAQA_WITH_CASADI_OCP_SUPPORTED} "ALPAQA_WITH_CASADI;ALPAQA_WITH_OCP" Off)
option(ALPAQA_WITH_DL
    "Enable dynamic loading of problems" On)
# Optimization and parallelization
option(ALPAQA_WITH_OPENMP
    "Enable parallelization using OpenMP" Off)
# Enable/disable floating point types
option(ALPAQA_WITH_QUAD_PRECISION
    "Support quad-precision floating point types" Off)
option(ALPAQA_WITH_SINGLE_PRECISION
    "Support single-precision floating point types" Off)
option(ALPAQA_WITH_LONG_DOUBLE
    "Support long double floating point types" Off)
# Developer options
option(ALPAQA_WARNINGS_AS_ERRORS
    "Enable -Werror or /WX" Off)
option(ALPAQA_WITH_COVERAGE
    "Generate coverage information" Off)
option(ALPAQA_ENABLE_PCH
    "Enable precompiled headers" Off)
set(ALPAQA_DOXYFILE "Doxyfile" CACHE FILEPATH
    "The Doxyfile to use for the docs target")
option(ALPAQA_DEBUG_CHECKS_EIGEN
    "Initialize Eigen matrices to NaN, and raise an assertion error if memory (re)allocation is used in the inner loops of the algorithms" Off)
option(ALPAQA_DONT_PARALLELIZE_EIGEN
    "Add the EIGEN_DONT_PARALLELIZE option" On)
option(ALPAQA_WITH_BLAS
    "Enable use of BLAS and enable EIGEN_USE_BLAS" Off)
option(ALPAQA_NO_DLCLOSE
    "Never call dlclose when a shared library is no longer needed, leave it open instead. Useful for profiling using gperftools." Off)
option(ALPAQA_WITH_ACCURATE_BUILD_TIME
    "Update the build time on every build" On)

# Installation paths
include(GNUInstallDirs)
set(ALPAQA_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}"
    CACHE PATH "Installation directory for archives and libraries")
set(ALPAQA_INSTALL_CMAKEDIR "${ALPAQA_INSTALL_LIBDIR}/cmake/alpaqa"
    CACHE PATH "Installation directory for CMake configuration files")
set(ALPAQA_INSTALL_BINDIR "${CMAKE_INSTALL_BINDIR}"
    CACHE PATH "Installation directory for binaries and DLLs")
set(ALPAQA_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}"
    CACHE PATH "Installation directory for headers")
set(ALPAQA_INSTALL_DATADIR "${CMAKE_INSTALL_DATADIR}"
    CACHE PATH "Installation directory for read-only architecture-independent data (share)")
set(ALPAQA_INSTALL_ZSHCOMPLETEDIR "${ALPAQA_INSTALL_DATADIR}/zsh/site-functions"
    CACHE PATH "Installation directory for Zsh completion scripts")
set(ALPAQA_INSTALL_BASHCOMPLETEDIR "${ALPAQA_INSTALL_DATADIR}/bash-completion/completions"
    CACHE PATH "Installation directory for Bash completion scripts")
set(ALPAQA_INSTALL_PYSTUBSDIR ""
    CACHE PATH "Installation directory for Python stubs (relative to Python package)")
option(ALPAQA_STANDALONE
    "Install with relative RPATH to locate its own shared libraries" On)

# Compiler warnings
include(cmake/Warnings.cmake)
add_warnings_target(warnings ${ALPAQA_WARNINGS_AS_ERRORS})
add_library(alpaqa::warnings ALIAS warnings)

# Compiler options
if (NOT DEFINED CMAKE_CXX_EXTENSIONS)
    set(CMAKE_CXX_EXTENSIONS ${ALPAQA_WITH_QUAD_PRECISION})
endif()
if (NOT DEFINED CMAKE_C_EXTENSIONS)
    set(CMAKE_C_EXTENSIONS ${ALPAQA_WITH_QUAD_PRECISION})
endif()
set(CMAKE_POSITION_INDEPENDENT_CODE On) # For Examples, Python/Matlab interfaces
if (MSVC)
    add_compile_options("/bigobj")
endif()
set(CMAKE_RELEASE_POSTFIX "")
set(CMAKE_DEBUG_POSTFIX "_d")
set(CMAKE_RELWITHDEBINFO_POSTFIX "_rd")
set(CMAKE_MINSIZEREL_POSTFIX "_rs")

# Sanitizers
include(cmake/Sanitizers.cmake)

# Coverage
if (ALPAQA_WITH_COVERAGE)
    add_custom_target(coverage
        ${CMAKE_CURRENT_LIST_DIR}/scripts/coverage.sh
        ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}
        ${CMAKE_CXX_COMPILER}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        USES_TERMINAL)
    add_compile_options("--coverage")
    add_link_options("--coverage")
    add_dependencies(coverage alpaqa::tests)
endif()

# Documentation
find_package(Doxygen)
if (DOXYGEN_FOUND)
    add_custom_target(docs Doxygen::doxygen ${ALPAQA_DOXYFILE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/doxygen)
endif()

# Save the project version to a generated file
configure_file(cmake/alpaqa-version.h.in
    ${PROJECT_BINARY_DIR}/include/alpaqa-version.h @ONLY)
include_directories(${PROJECT_BINARY_DIR}/include)

# Libraries
add_subdirectory(src)

# Interfaces to other languages
add_subdirectory(interfaces)

# Tests
if (ALPAQA_WITH_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

# Examples
if (ALPAQA_WITH_EXAMPLES)
    add_subdirectory(examples)
endif()

# Packaging
include(cmake/Packaging.cmake)
