cmake_minimum_required(VERSION 2.6)
include_directories(${QPOASES_INCLUDE_DIR})

set(QPOASES_INTERFACE_SRCS
  qpoases_internal.hpp
  qpoases_internal.cpp
  qpoases_solver.hpp
  qpoases_solver.cpp)

casadi_library(casadi_qpoases_interface ${QPOASES_INTERFACE_SRCS})
target_link_libraries(casadi_qpoases_interface ${QPOASES_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/qpoases
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_qpoases_interface PROPERTIES COMPILE_FLAGS "-Wno-overloaded-virtual")
endif()
