cmake_minimum_required(VERSION 2.0)

include_directories(../)
find_package(GSL REQUIRED QUIET)
include_directories(${GSL_INCLUDE_DIRS})

set(GSL_INTERFACE_SRCS
  gsl_internal.hpp
  gsl_internal.cpp
  gsl_integrator.hpp
  gsl_integrator.cpp
)

add_library(gsl_interface STATIC ${GSL_INTERFACE_SRCS})
if(ENABLE_SHARED)
add_library(gsl_interface_shared SHARED ${GSL_INTERFACE_SRCS})
endif(ENABLE_SHARED)

