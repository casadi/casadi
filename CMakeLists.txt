cmake_minimum_required(VERSION 3.16)

project(PANOC-ALM VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# find_package(GTest REQUIRED PATHS "$ENV{HOME}/.local/lib/cmake" NO_DEFAULT_PATH)
find_package(GTest REQUIRED NO_CMAKE_SYSTEM_PATH)
include(GoogleTest)
enable_testing()

# Set the compiler warnings
set(CMAKE_CPP_FLAGS "${CMAKE_CPP_FLAGS} \
    -Wmissing-include-dirs")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -Wall -Wextra -Werror -pedantic -pedantic-errors -fdiagnostics-show-option \
    -Wno-unknown-pragmas")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -Wall -Wextra -Werror -pedantic -pedantic-errors -fdiagnostics-show-option")
    
# Build-type specific flags
set(CMAKE_CXX_FLAGS_DEBUG "-ggdb3 -O0 -D_GLIBCXX_DEBUG")
set(CMAKE_C_FLAGS_DEBUG   "-ggdb3 -O0 -D_GLIBCXX_DEBUG")
set(CMAKE_DEBUG_POSTFIX "-debug")

set(CMAKE_CXX_FLAGS_ASAN "-g3 -O0 -D_GLIBCXX_DEBUG \
    -fsanitize=address,leak,undefined,pointer-compare,pointer-subtract")
set(CMAKE_C_FLAGS_ASAN   "-g3 -O0 -D_GLIBCXX_DEBUG \
    -fsanitize=address,leak,undefined,pointer-compare,pointer-subtract")
set(CMAKE_ASAN_POSTFIX "-asan")

set(CMAKE_CXX_FLAGS_TSAN "-g3 -O0 -D_GLIBCXX_DEBUG \
    -fsanitize=thread,undefined")
set(CMAKE_C_FLAGS_TSAN   "-g3 -O0 -D_GLIBCXX_DEBUG \
    -fsanitize=thread,undefined")
set(CMAKE_TSAN_POSTFIX "-tsan")

# Coverage
set(CMAKE_C_FLAGS_COVERAGE "${CMAKE_C_FLAGS_ASAN} \
    --coverage -fno-inline")
set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_ASAN} \
    --coverage -fno-inline")
set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "--coverage")
set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)

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
