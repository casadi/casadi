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
  ipopt_qp_internal.hpp
  ipopt_qp_internal.cpp
  ipopt_qp_solver.hpp
  ipopt_qp_solver.cpp
)


if(ENABLE_STATIC)
add_library(ipopt_interface STATIC ${IPOPT_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(ipopt_interface SHARED ${IPOPT_INTERFACE_SRCS})
endif(ENABLE_SHARED)
