cmake_minimum_required(VERSION 2.0)
find_package(ACADO REQUIRED)
include_directories(${ACADO_INCLUDE_DIRS})

set(ACADO_INTERFACE_SRCS
  acado_function.hpp
  acado_function.cpp
  acado_interface.hpp
  acado_interface.cpp
  acado_internal.hpp
  acado_internal.cpp
  acado_integrator_backend.hpp
  acado_integrator_backend.cpp
)

add_library(acado_interface STATIC ${ACADO_INTERFACE_SRCS})
if(ENABLE_SHARED)
add_library(acado_interface_shared SHARED ${ACADO_INTERFACE_SRCS})
endif(ENABLE_SHARED)

