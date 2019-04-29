cmake_minimum_required(VERSION 2.8.6)

include_directories(../)

casadi_plugin(Conic qpalm
  qpalm_interface.hpp
  qpalm_interface.cpp
  qpalm_interface_meta.cpp)

casadi_plugin_link_libraries(Conic qpalm qpalm)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_conic_osqp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-unknown-warning-option")
endif()
