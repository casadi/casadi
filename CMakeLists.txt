cmake_minimum_required(VERSION 3.24)
project(alpaqa VERSION 1.0.0 LANGUAGES CXX)
set(PY_VERSION_SUFFIX "a7")

# Look for a C and Fortran compiler (optional)
include(CheckLanguage)
set(ALPAQA_HAVE_FORTRAN Off)
check_language(Fortran)
if (CMAKE_Fortran_COMPILER)
    set(ALPAQA_HAVE_FORTRAN On)
    enable_language(Fortran)
endif()
set(ALPAQA_HAVE_C Off)
check_language(C)
if (CMAKE_C_COMPILER)
    set(ALPAQA_HAVE_C On)
    enable_language(C)
endif()

# Options
option(ALPAQA_WITH_QUAD_PRECISION
    "Support quad-precision floating point types" Off)
option(ALPAQA_WITH_PYTHON
    "Build the Python bindings" Off)
option(ALPAQA_WARNINGS_AS_ERRORS
    "Enable -Werror or /WX" Off)
option(ALPAQA_WITH_TESTS
    "Build the tests" On)
option(ALPAQA_WITH_EXAMPLES
    "Build the examples" On)
option(ALPAQA_WITH_DRIVERS
    "Build the solver drivers" On)
option(ALPAQA_WITH_COVERAGE
    "Generate coverage information" Off)
option(ALPAQA_WITH_CASADI
    "Build the CasADi loader" On)
set(ALPAQA_DOXYFILE "Doxyfile" CACHE STRING
    "The Doxyfile to use for the docs target")
option(ALPAQA_ENABLE_PCH
    "Enable precompiled headers" On)

# Compiler warnings
include(cmake/Warnings.cmake)
add_warnings_target(warnings ${ALPAQA_WARNINGS_AS_ERRORS})
add_library(alpaqa::warnings ALIAS warnings)

# Compiler options
if (ALPAQA_WITH_PYTHON)
    set(CMAKE_POSITION_INDEPENDENT_CODE On)
endif()
if (MSVC)
    set(CMAKE_CXX_FLAGS "/utf-8 /Dnot=! /D__PRETTY_FUNCTION__=__FUNCSIG__ ${CMAKE_CXX_FLAGS} /bigobj /permissive-")
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
        ${CMAKE_CXX_COMPILER_ID}
        ${CMAKE_CXX_COMPILER_VERSION}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage -fno-inline")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -fno-inline")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
    add_dependencies(coverage alpaqa::tests)
endif()

# Documentation
find_package(Doxygen)
if (DOXYGEN_FOUND)
    add_custom_target(docs
        Doxygen::doxygen ${ALPAQA_DOXYFILE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/doxygen)
endif()

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
include(InstallRequiredSystemLibraries)
set(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_PATCH}${PY_VERSION_SUFFIX}")
set(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/README.md")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Pieter Pas <pieter.p.dev@outlook.com>")
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "https://github.com/kul-optec/alpaqa")
# Package components
set(CPACK_COMPONENTS_ALL "shlib;dev;dl_dev;debug;bin")
set(CPACK_DEB_COMPONENT_INSTALL On)
set(CPACK_DEBIAN_SHLIB_PACKAGE_NAME "libalpaqa")
set(CPACK_DEBIAN_DEV_PACKAGE_NAME "libalpaqa-dev")
set(CPACK_DEBIAN_DL_DEV_PACKAGE_NAME "libalpaqa-dl-dev")
set(CPACK_DEBIAN_DEBUG_PACKAGE_NAME "libalpaqa-debug")
set(CPACK_DEBIAN_BIN_PACKAGE_NAME "libalpaqa-tools")
# Package dependencies
set(CPACK_DEBIAN_SHLIB_PACKAGE_DEPENDS "libc6 (>= 2.17)")
set(CPACK_DEBIAN_DEV_PACKAGE_DEPENDS "libalpaqa (= ${PROJECT_VERSION}${PY_VERSION_SUFFIX}), libalpaqa-dl-dev (= ${PROJECT_VERSION}${PY_VERSION_SUFFIX})")
set(CPACK_DEBIAN_DEV_PACKAGE_SUGGESTS "libeigen3-dev (>= 3.4.0), libeigen3-dev (<< 3.4.1)")
set(CPACK_DEBIAN_DL_DEV_PACKAGE_DEPENDS "")
set(CPACK_DEBIAN_DL_DEV_PACKAGE_SUGGESTS "libalpaqa-dev (= ${PROJECT_VERSION}${PY_VERSION_SUFFIX})")
set(CPACK_DEBIAN_DEBUG_PACKAGE_DEPENDS "libalpaqa (= ${PROJECT_VERSION}${PY_VERSION_SUFFIX})")
set(CPACK_DEBIAN_BIN_PACKAGE_DEPENDS "libalpaqa (= ${PROJECT_VERSION}${PY_VERSION_SUFFIX})")
include(CPack)
