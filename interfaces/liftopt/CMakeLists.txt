cmake_minimum_required(VERSION 2.0)
find_package(LIFTOPT REQUIRED QUIET)
include_directories(${LIFTOPT_INCLUDE_DIRS})

SET(LIFTOPT_INTERFACE_SRCS
  liftopt_solver.hpp
  liftopt_solver.cpp
  liftopt_internal.hpp
  liftopt_internal.cpp
  casadi_lifter.hpp
  casadi_lifter.cpp
)

add_library(liftopt_interface STATIC ${LIFTOPT_INTERFACE_SRCS})
if(ENABLE_SHARED)
add_library(liftopt_interface_shared SHARED ${LIFTOPT_INTERFACE_SRCS})
endif(ENABLE_SHARED)
