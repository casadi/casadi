cmake_minimum_required(VERSION 2.8.6)

set(CMAKE_CXX_STANDARD 17)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(../)

casadi_plugin(Conic piqp
  piqp_interface.hpp
  piqp_interface.cpp
  piqp_interface_meta.cpp)

casadi_plugin_link_libraries(Conic piqp piqpcstatic Eigen3::Eigen)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_conic_piqp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-unknown-warning-option")
endif()

set_property(TARGET casadi_conic_piqp PROPERTY PKG_CONFIG_REQUIRES piqp)
