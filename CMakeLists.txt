cmake_minimum_required(VERSION 3.17)

project(alpaqa VERSION 0.0.1 LANGUAGES C CXX)

if (DEFINED VERIFY_VERSION)
    if (NOT "${VERIFY_VERSION}" MATCHES "^${PROJECT_VERSION}.*")
        message(FATAL_ERROR "Version number does not match "
                             "(${VERIFY_VERSION} - ${PROJECT_VERSION}).")
    endif()
endif()

if (NOT WIN32)
    include(CheckLanguage)
    check_language(Fortran)
    if (CMAKE_Fortran_COMPILER)
        enable_language(Fortran)
    endif()
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

if (DEFINED ENV{VIRTUAL_ENV})
    list(APPEND CMAKE_PREFIX_PATH "$ENV{VIRTUAL_ENV}")
elseif (EXISTS "${CMAKE_CURRENT_LIST_DIR}/py-venv")
    list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_LIST_DIR}/py-venv")
endif()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# Import file with some custom compiler warnings and debug/sanitizer flags 
include(cmake/GlobalCompilerOptions.cmake)
# Macro for generating objective/constraint functions and derivatives in Python
# using CasADi
include(cmake/CasADiFunctionCodeGenPython.cmake)
# Macro for generating CUTEst test problems
include(cmake/CUTEstSIFProblem.cmake)

# ----

# Google Test
find_package(GTest REQUIRED NO_CMAKE_SYSTEM_PATH)
include(GoogleTest)
enable_testing()

# Eigen
find_package(Eigen3 REQUIRED)

# CasADi
find_package(casadi)

# CUTEst
find_package(CUTEst)

# LBFGSpp
find_package(LBFGSpp)

# ----

add_subdirectory(src)
if (NOT SKBUILD)
    add_subdirectory(test)
    add_subdirectory(examples)
endif()

# ----

# Add documentation targets
add_custom_target(documentation-doxygen
    doxygen
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/doxygen)
add_custom_target(documentation-doxygen-breathe
    doxygen Doxyfile.breathe
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/doxygen)
add_custom_target(documentation-sphinx-doctest
    sphinx-build -b doctest
        ${CMAKE_CURRENT_LIST_DIR}/sphinx/source 
        ${CMAKE_CURRENT_LIST_DIR}/docs/Sphinx
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_dependencies(documentation-sphinx-doctest documentation-doxygen-breathe)
add_custom_target(documentation-sphinx
    sphinx-build -M html
        ${CMAKE_CURRENT_LIST_DIR}/sphinx/source 
        ${CMAKE_CURRENT_LIST_DIR}/docs/Sphinx
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_dependencies(documentation-sphinx documentation-doxygen-breathe)
add_custom_target(documentation)
add_dependencies(documentation documentation-doxygen documentation-sphinx)

# Add coverage target
if (CMAKE_BUILD_TYPE STREQUAL "Coverage")
    add_custom_target(coverage 
        ${CMAKE_CURRENT_LIST_DIR}/scripts/coverage.sh
        ${CMAKE_CXX_COMPILER_ID}
        ${CMAKE_CXX_COMPILER_VERSION}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

# ----

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${alpaqa_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${alpaqa_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${alpaqa_VERSION_PATCH}")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
SET(CPACK_GENERATOR "DEB;TGZ")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Pieter P")
include(CPack)
