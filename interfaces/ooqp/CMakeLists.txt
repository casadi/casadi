cmake_minimum_required(VERSION 2.0)

include_directories(../)
find_package(OOQP REQUIRED QUIET)
include_directories(${OOQP_INCLUDE_DIRS})

set(OOQP_INTERFACE_SRCS
  ooqp_internal.hpp
  ooqp_internal.cpp
  ooqp_solver.hpp
  ooqp_solver.cpp
)

add_library(ooqp_interface STATIC ${OOQP_INTERFACE_SRCS})
if(ENABLE_SHARED)
add_library(ooqp_interface_shared SHARED ${OOQP_INTERFACE_SRCS})
endif(ENABLE_SHARED)

