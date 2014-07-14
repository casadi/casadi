cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${OOQP_INCLUDE_DIRS})

casadi_library(casadi_qpsolver_ooqp
  ooqp_internal.hpp
  ooqp_internal.cpp)

target_link_libraries(casadi_qpsolver_ooqp ${OOQP_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/ooqp
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)


