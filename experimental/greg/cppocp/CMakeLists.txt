cmake_minimum_required(VERSION 2.0)
project(cppocp)

SET(CPPOCP_SRCS
  Ddp.cpp   Ddp.hpp
  Lqr.cpp   Lqr.hpp
  Ocp.cpp   Ocp.hpp
  MultipleShooting.cpp   MultipleShooting.hpp
  Ode.cpp   Ode.hpp
)

add_library(cppocp STATIC ${CPPOCP_SRCS})
if(ENABLE_SHARED)
add_library(cppocp_shared SHARED ${CPPOCP_SRCS})
endif(ENABLE_SHARED)
