cmake_minimum_required(VERSION 3.10.2)

option(USE_SYSTEM_WISE_OSQP "Use a systemwise installation of osqp" OFF)

include_directories(../)

casadi_plugin(Conic osqp
  osqp_interface.hpp
  osqp_interface.cpp
  osqp_interface_meta.cpp)

# Add the possibility to use a systemwise installation of casadi
if(USE_SYSTEM_WISE_OSQP)
  find_package(osqp REQUIRED)
  casadi_plugin_link_libraries(Conic osqp osqp::osqp)
else()
  casadi_plugin_link_libraries(Conic osqp osqp)
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_conic_osqp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-unknown-warning-option")
endif()
