cmake_minimum_required(VERSION 2.0)
project(cppocp)

SET(CPPOCP_SRCS
  Ocp.cpp   Ocp.hpp
  OcpMultipleShooting.cpp   OcpMultipleShooting.hpp
  Ode.cpp   Ode.hpp
)

add_library(cppocp STATIC ${CPPOCP_SRCS})
if(ENABLE_SHARED)
add_library(cppocp_shared SHARED ${CPPOCP_SRCS})
endif(ENABLE_SHARED)
