cmake_minimum_required(VERSION 2.0)
find_package(IPOPT REQUIRED QUIET)
include_directories(${IPOPT_INCLUDE_DIR})

SET(IPOPT_INTERFACE_SRCS
  ipopt_solver.hpp
  ipopt_solver.cpp
  ipopt_internal.hpp
  ipopt_internal.cpp
  ipopt_nlp.hpp
  ipopt_nlp.cpp
#   ipopt_solver_c.h
#   ipopt_solver_c.cpp
)


add_library(ipopt_interface STATIC ${IPOPT_INTERFACE_SRCS})
if(ENABLE_SHARED)
add_library(ipopt_interface_shared SHARED ${IPOPT_INTERFACE_SRCS})
endif(ENABLE_SHARED)


add_subdirectory(test)
