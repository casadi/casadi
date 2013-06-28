cmake_minimum_required(VERSION 2.6)
include_directories(${DSDP_INCLUDE_DIR})

set(DSDP_INTERFACE_SRCS
  dsdp_internal.hpp
  dsdp_internal.cpp
  dsdp_solver.hpp
  dsdp_solver.cpp
)

if(ENABLE_STATIC)
add_library(casadi_dsdp_interface STATIC ${DSDP_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_dsdp_interface SHARED ${DSDP_INTERFACE_SRCS})
endif(ENABLE_SHARED)
install(TARGETS casadi_dsdp_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/interfaces/dsdp
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)
