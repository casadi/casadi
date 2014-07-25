cmake_minimum_required(VERSION 2.6)
include_directories(${KNITRO_INCLUDE_DIR})

casadi_plugin(NlpSolver knitro
  knitro_internal.hpp
  knitro_internal.cpp
  knitro_internal_meta.cpp)

casadi_plugin_link_libraries(NlpSolver knitro ${KNITRO_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/knitro
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
