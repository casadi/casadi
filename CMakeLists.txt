cmake_minimum_required(VERSION 3.16)

project(PANOC-ALM VERSION 0.0.1 LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_LIST_DIR}/py-venv")
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

# ----

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(examples)

# ----

# Add export documentation target
add_custom_target(documentation 
    doxygen
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/doxygen)

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
set(CPACK_PACKAGE_VERSION_MAJOR "${PANOC-ALM_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PANOC-ALM_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${PANOC-ALM_VERSION_PATCH}")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
SET(CPACK_GENERATOR "DEB;TGZ")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Pieter P")
include(CPack)
