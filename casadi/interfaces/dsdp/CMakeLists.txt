cmake_minimum_required(VERSION 2.6)
include_directories(${DSDP_INCLUDE_DIR})

casadi_library(casadi_sdpsolver_dsdp
  dsdp_internal.hpp
  dsdp_internal.cpp)

target_link_libraries(casadi_sdpsolver_dsdp ${DSDP_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/dsdp
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
