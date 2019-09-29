cmake_minimum_required(VERSION 2.8.6)

# Get package info using pkg-config
find_package(PkgConfig)
pkg_search_module(CASADI casadi)

if(NOT CASADI_LIBRARIES)
  message(SEND_ERROR "Could not find casadi.pc - did you set PKG_CONFIG_PATH?")
endif()

message("CASADI_LIBRARY_DIRS: ${CASADI_LIBRARY_DIRS}")
message("CASADI_LIBRARIES: ${CASADI_LIBRARIES}")
message("CASADI_INCLUDEDIR: ${CASADI_INCLUDEDIR}")

link_directories(${CASADI_LIBRARY_DIRS})
add_executable(casadi_demo casadi_demo.cpp)


target_link_libraries(casadi_demo ${CASADI_LIBRARIES})
target_include_directories(casadi_demo PRIVATE ${CASADI_INCLUDEDIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
