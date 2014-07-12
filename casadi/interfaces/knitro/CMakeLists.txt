cmake_minimum_required(VERSION 2.6)
include_directories(${KNITRO_INCLUDE_DIR})

casadi_library(casadi_nlpsolver_knitro
  knitro_internal.hpp
  knitro_internal.cpp)

target_link_libraries(casadi_nlpsolver_knitro ${KNITRO_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/knitro
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
