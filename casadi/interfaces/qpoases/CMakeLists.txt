cmake_minimum_required(VERSION 2.6)
include_directories(${QPOASES_INCLUDE_DIR})

casadi_library(casadi_qpsolver_qpoases
  qpoases_internal.hpp
  qpoases_internal.cpp)

target_link_libraries(casadi_qpsolver_qpoases ${QPOASES_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/qpoases
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_qpsolver_qpoases PROPERTIES COMPILE_FLAGS "-Wno-overloaded-virtual")
endif()
