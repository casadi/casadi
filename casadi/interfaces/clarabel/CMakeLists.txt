cmake_minimum_required(VERSION 3.10.2)


casadi_plugin(Conic clarabel
  clarabel_interface.hpp
  clarabel_interface.cpp
  clarabel_interface_meta.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/clarabel_runtime_str.h"
  )

casadi_plugin_link_libraries(Conic clarabel clarabel)

# Add the runtime sources to internal (installed separately)
set(RUNTIME_SRC "clarabel_runtime.hpp")

string (REPLACE ";" "$<SEMICOLON>" ESCAPED_RUNTIME_SRC "${RUNTIME_SRC}")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/clarabel_runtime_str.h"
  COMMAND  ${CMAKE_COMMAND} -D OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/clarabel_runtime_str.h" -D SOURCES="${RUNTIME_SRC}" -P "${CMAKE_CURRENT_SOURCE_DIR}/../../generate_runtime.cmake"
  DEPENDS ${RUNTIME_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_include_directories(casadi_conic_clarabel
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)