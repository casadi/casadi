cmake_minimum_required(VERSION 2.8.6)

set(CMAKE_CXX_STANDARD 17)

casadi_plugin(Conic proxqp
  proxqp_interface.hpp
  proxqp_interface.cpp
  proxqp_interface_meta.cpp)


casadi_plugin_link_libraries(Conic proxqp)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set_target_properties(casadi_conic_proxqp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-unknown-warning-option")
endif()

target_include_directories(casadi_conic_proxqp PUBLIC ${proxsuite_INCLUDE_DIRS})
target_link_libraries(casadi_conic_proxqp proxsuite)
