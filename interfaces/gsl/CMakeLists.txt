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

if(ENABLE_STATIC)
add_library(gsl_interface STATIC ${GSL_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(gsl_interface SHARED ${GSL_INTERFACE_SRCS})
endif(ENABLE_SHARED)

