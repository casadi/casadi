cmake_minimum_required(VERSION 2.8.6)
include_directories(${QPOASES_INCLUDE_DIR})

casadi_plugin(QpSolver qp qpoases
  qpoases_interface.hpp
  qpoases_interface.cpp
  qpoases_interface_meta.cpp)

casadi_plugin_link_libraries(qp qpoases ${QPOASES_LIBRARIES})

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_qp_qpoases PROPERTIES COMPILE_FLAGS "-Wno-overloaded-virtual")
endif()
