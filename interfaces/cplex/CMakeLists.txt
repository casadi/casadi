cmake_minimum_required(VERSION 2.0)
find_package(CPLEX REQUIRED QUIET)
include_directories(${CPLEX_INCLUDE_DIR})

SET(CPLEX_INTERFACE_SRCS
  cplex_solver.hpp
  cplex_solver.cpp
  cplex_internal.hpp
  cplex_internal.cpp
)

if(ENABLE_STATIC)
add_library(cplex_interface STATIC ${CPLEX_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(cplex_interface SHARED ${CPLEX_INTERFACE_SRCS})
endif(ENABLE_SHARED)
