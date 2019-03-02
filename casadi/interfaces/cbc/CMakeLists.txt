cmake_minimum_required(VERSION 2.8.6)

include_directories(${CBC_INCLUDE_DIRS})
link_directories(${CBC_LIBRARY_DIRS})

casadi_plugin(Conic cbc
  cbc_interface.hpp
  cbc_interface.cpp
  cbc_interface_meta.cpp)
casadi_plugin_link_libraries(Conic cbc ${CBC_LIBRARIES})

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_conic_cbc PROPERTIES COMPILE_FLAGS "-Wno-misleading-indentation -Wno-unknown-warning-option")
endif()
