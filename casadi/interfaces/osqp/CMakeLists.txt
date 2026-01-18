cmake_minimum_required(VERSION 3.10.2)

include_directories(../)

casadi_plugin(Conic osqp
  osqp_interface.hpp
  osqp_interface.cpp
  osqp_interface_meta.cpp)

casadi_plugin_link_libraries(Conic osqp osqp::osqp)


if(NOT DEFINED WITH_OSQP_V1)
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/try-osqp-v1.cpp"
"#include <osqp_api_functions.h>

int main()
{
  // This function is only available in OSQP >= 1.0.0
  OSQPCscMatrix_set_data(nullptr, 0, 0, 0, nullptr, 0, 0);
  return 0;
}")

  try_compile(WITH_OSQP_V1 ${CMAKE_CURRENT_BINARY_DIR} "${CMAKE_CURRENT_BINARY_DIR}/try-osqp-v1.cpp" LINK_LIBRARIES osqp::osqp)
endif()

if(WITH_OSQP_V1)
  target_compile_definitions(casadi_conic_osqp PRIVATE WITH_OSQP_V1)
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_conic_osqp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-unknown-warning-option")
endif()
