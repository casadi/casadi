cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${OOQP_INCLUDE_DIRS})

set(OOQP_INTERFACE_SRCS
  ooqp_internal.hpp
  ooqp_internal.cpp
  ooqp_solver.hpp
  ooqp_solver.cpp)

casadi_library(casadi_ooqp_interface ${OOQP_INTERFACE_SRCS})
target_link_libraries(casadi_ooqp_interface ${OOQP_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/ooqp
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)


