cmake_minimum_required(VERSION 2.8.12)

include_directories(../)

casadi_plugin(Conic qpswift
  qpswift_interface.hpp
  qpswift_interface.cpp
  qpswift_interface_meta.cpp)

casadi_plugin_link_libraries(Conic qpswift qpswift)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_conic_qpswift PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-unknown-warning-option")
endif()
