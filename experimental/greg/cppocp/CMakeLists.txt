cmake_minimum_required(VERSION 2.6)
project(cppocp)

set(CPPOCP_SRCS
  Ddp.cpp   Ddp.hpp
  Lqr.cpp   Lqr.hpp
  Ocp.cpp   Ocp.hpp
  MultipleShooting.cpp   MultipleShooting.hpp
  Ode.cpp   Ode.hpp)

if(ENABLE_STATIC)
add_library(cppocp STATIC ${CPPOCP_SRCS})
endif()
if(ENABLE_SHARED)
add_library(cppocp SHARED ${CPPOCP_SRCS})
endif()
