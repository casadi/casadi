cmake_minimum_required(VERSION 2.0)
find_package(KNITRO REQUIRED QUIET)
include_directories(${KNITRO_INCLUDE_DIR})

SET(KNITRO_INTERFACE_SRCS
  knitro_solver.hpp
  knitro_solver.cpp
  knitro_internal.hpp
  knitro_internal.cpp
)


add_library(knitro_interface STATIC ${KNITRO_INTERFACE_SRCS})
if(ENABLE_SHARED)
add_library(knitro_interface_shared SHARED ${KNITRO_INTERFACE_SRCS})
endif(ENABLE_SHARED)
