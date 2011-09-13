cmake_minimum_required(VERSION 2.0)
find_package(ACADO REQUIRED)
include_directories(${ACADO_INCLUDE_DIRS})

set(ACADO_INTERFACE_SRCS
  acado_function.hpp
  acado_function.cpp
  acado_ocp_interface.hpp
  acado_ocp_interface.cpp
  acado_ocp_internal.hpp
  acado_ocp_internal.cpp
  acado_integrator_backend.hpp
  acado_integrator_backend.cpp
)

if(ENABLE_STATIC)
add_library(acado_interface STATIC ${ACADO_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(acado_interface SHARED ${ACADO_INTERFACE_SRCS})
endif(ENABLE_SHARED)
install(TARGETS acado_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./
   DESTINATION include/casadi/interfaces/acado
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)

