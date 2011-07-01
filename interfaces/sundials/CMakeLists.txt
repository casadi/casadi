cmake_minimum_required(VERSION 2.0)

include_directories(../)
find_package(SUNDIALS REQUIRED QUIET)
include_directories(${SUNDIALS_INCLUDE_DIR})
if(PRECOND_TEST)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPRECOND_TEST")
endif(PRECOND_TEST)

set(SUNDIALS_INTERFACE_SRCS
  cvodes_internal.hpp
  cvodes_internal.cpp
  cvodes_integrator.hpp
  cvodes_integrator.cpp
  idas_internal.hpp
  idas_internal.cpp
  idas_integrator.hpp
  idas_integrator.cpp
  kinsol_internal.hpp
  kinsol_internal.cpp
  kinsol_solver.hpp
  kinsol_solver.cpp
)

if(ENABLE_STATIC)
add_library(sundials_interface STATIC ${SUNDIALS_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(sundials_interface SHARED ${SUNDIALS_INTERFACE_SRCS})
endif(ENABLE_SHARED)

