cmake_minimum_required(VERSION 3.10.2)

set(CMAKE_CXX_STANDARD 17)

casadi_plugin(Conic proxqp
  proxqp_interface.hpp
  proxqp_interface.cpp
  proxqp_interface_meta.cpp)


casadi_plugin_link_libraries(Conic proxqp proxqp)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set_target_properties(casadi_conic_proxqp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-unknown-warning-option -Wno-type-limits")
endif()