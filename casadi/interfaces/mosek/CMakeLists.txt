cmake_minimum_required(VERSION 3.10.2)

include_directories(../)

casadi_plugin(Conic mosek
  mosek_interface.hpp
  mosek_interface.cpp
  mosek_interface_meta.cpp
  mosek_runtime.hpp
  "${CMAKE_CURRENT_BINARY_DIR}/mosek_runtime_str.h"
  )

casadi_plugin_link_libraries(Conic mosek mosek::mosek)

# Generate runtime string header for codegen
set(RUNTIME_SRC "mosek_runtime.hpp")

string (REPLACE ";" "$<SEMICOLON>" ESCAPED_RUNTIME_SRC "${RUNTIME_SRC}")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/mosek_runtime_str.h"
  COMMAND  ${CMAKE_COMMAND} -D OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/mosek_runtime_str.h" -D SOURCES="${RUNTIME_SRC}" -P "${CMAKE_CURRENT_SOURCE_DIR}/../../generate_runtime.cmake"
  DEPENDS ${RUNTIME_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_include_directories(casadi_conic_mosek
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
