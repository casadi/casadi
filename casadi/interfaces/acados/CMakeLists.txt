cmake_minimum_required(VERSION 3.16.3)

include_directories(../)

casadi_plugin(Integrator acados
  acados_integrator.hpp
  acados_integrator.cpp
  acados_integrator_meta.cpp
  acados_chain.hpp
  acados_sim.hpp
  acados_runtime.hpp
  "${CMAKE_CURRENT_BINARY_DIR}/acados_runtime_str.h")

# Runtime embedded in generated code
set(RUNTIME_SRC "acados_sim.hpp" "acados_chain.hpp" "acados_runtime.hpp")

string (REPLACE ";" "$<SEMICOLON>" ESCAPED_RUNTIME_SRC "${RUNTIME_SRC}")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/acados_runtime_str.h"
  COMMAND  ${CMAKE_COMMAND} -D OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/acados_runtime_str.h" -D SOURCES="${ESCAPED_RUNTIME_SRC}" -P "${CMAKE_CURRENT_SOURCE_DIR}/../../generate_runtime.cmake"
  DEPENDS ${RUNTIME_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_include_directories(casadi_integrator_acados
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

casadi_plugin_link_libraries(Integrator acados acados)
